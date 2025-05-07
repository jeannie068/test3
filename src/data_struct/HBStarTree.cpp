/**
 * HBStarTree.cpp
 * 
 * Implementation of the HBStarTree class for analog placement with symmetry constraints.
 * The HB*-tree is a hierarchical framework that can simultaneously optimize the placement
 * with both symmetry islands and non-symmetric modules.
 */

#include "HBStarTree.hpp"
#include <algorithm>
#include <queue>
#include <iostream>
#include <limits>
#include <stack>
using namespace std;

/**
 * Constructor
 */
HBStarTree::HBStarTree()
    : root(nullptr),
      horizontalContour(make_shared<Contour>()),
      verticalContour(make_shared<Contour>()),
      totalArea(0),
      isPacked(false) {
}

HBStarTree::~HBStarTree() {

}

/* Add a module to the tree */
void HBStarTree::addModule(shared_ptr<Module> module) {
    if (!module) return;
    
    // Add the module to module map
    modules[module->getName()] = module;
}

/* Add a symmetry group to the tree */
void HBStarTree::addSymmetryGroup(shared_ptr<SymmetryGroup> group) {
    if (!group) return;
    
    // Add the symmetry group to list
    symmetryGroups.push_back(group);
}

/* Constructs the symmetry islands for each symmetry group */
void HBStarTree::constructSymmetryIslands() {
    // Create an ASF-B*-tree for each symmetry group
    for (const auto& group : symmetryGroups) {
        auto asfTree = make_shared<ASFBStarTree>(group);
        
        // Add modules that belong to this symmetry group
        for (const auto& pair : group->getSymmetryPairs()) {
            if (modules.find(pair.first) != modules.end()) {
                asfTree->addModule(modules[pair.first]);
            }
            if (modules.find(pair.second) != modules.end()) {
                asfTree->addModule(modules[pair.second]);
            }
        }
        
        for (const auto& moduleName : group->getSelfSymmetric()) {
            if (modules.find(moduleName) != modules.end()) {
                asfTree->addModule(modules[moduleName]);
            }
        }
        
        // Construct the initial ASF-B*-tree
        asfTree->constructInitialTree();
        
        // Create a hierarchy node for this symmetry group
        auto hierarchyNode = make_shared<HBStarTreeNode>(HBNodeType::HIERARCHY, group->getName());
        hierarchyNode->setASFTree(asfTree);
        
        // Add the hierarchy node to our map
        symmetryGroupNodes[group->getName()] = hierarchyNode;
    }
}

/* Constructs the initial tree structure */
void HBStarTree::constructInitialTreeStructure() {
    // Collect all non-symmetry modules
    vector<string> nonSymmetryModules;
    
    // Create a set of all modules in symmetry groups
    set<string> symmetryModules;
    for (const auto& group : symmetryGroups) {
        for (const auto& pair : group->getSymmetryPairs()) {
            symmetryModules.insert(pair.first);
            symmetryModules.insert(pair.second);
        }
        for (const auto& moduleName : group->getSelfSymmetric()) {
            symmetryModules.insert(moduleName);
        }
    }
    
    // Find non-symmetry modules
    for (const auto& pair : modules) {
        if (symmetryModules.find(pair.first) == symmetryModules.end()) {
            nonSymmetryModules.push_back(pair.first);
        }
    }
    
    // Sort non-symmetry modules by area (largest first) for better initial placement
    sort(nonSymmetryModules.begin(), nonSymmetryModules.end(), 
              [this](const string& a, const string& b) {
                  return modules[a]->getArea() > modules[b]->getArea();
              });
    
    // Create nodes for non-symmetry modules
    for (const auto& moduleName : nonSymmetryModules) {
        auto node = make_shared<HBStarTreeNode>(HBNodeType::MODULE, moduleName);
        moduleNodes[moduleName] = node;
    }
    
    // Create the initial tree - simplest approach is a left-skewed tree
    if (!symmetryGroupNodes.empty() || !moduleNodes.empty()) {
        // Start with the first symmetry group as root, if any
        if (!symmetryGroupNodes.empty()) {
            root = symmetryGroupNodes.begin()->second;
            auto current = root;
            
            // Add remaining symmetry groups
            for (auto it = next(symmetryGroupNodes.begin()); it != symmetryGroupNodes.end(); ++it) {
                current->setLeftChild(it->second);
                it->second->setParent(current);
                current = it->second;
            }
            
            // Add non-symmetry modules
            for (const auto& pair : moduleNodes) {
                current->setLeftChild(pair.second);
                pair.second->setParent(current);
                current = pair.second;
            }
        }
        // Otherwise, start with the first non-symmetry module
        else if (!moduleNodes.empty()) {
            auto it = moduleNodes.begin();
            root = it->second;
            auto current = root;
            
            // Add remaining non-symmetry modules
            for (++it; it != moduleNodes.end(); ++it) {
                current->setLeftChild(it->second);
                it->second->setParent(current);
                current = it->second;
            }
        }
    }
}

/* Clears the tree */
void HBStarTree::clearTree() {
    root = nullptr;
    moduleNodes.clear();
    symmetryGroupNodes.clear();
    isPacked = false;
}

/* Constructs an initial HB*-tree */
void HBStarTree::constructInitialTree() {
    // Clear any existing tree
    clearTree();
    
    // First, construct symmetry islands for each symmetry group
    constructSymmetryIslands();
    
    // Then, construct the initial tree structure
    constructInitialTreeStructure();
}

/**
 * Finds the nearest contour node for a given node
 */
shared_ptr<HBStarTreeNode> HBStarTree::findNearestContourNode(shared_ptr<HBStarTreeNode> node) const {
    if (!node || !root) return nullptr;
    
    // Use BFS to find the nearest contour node
    queue<shared_ptr<HBStarTreeNode>> queue;
    queue.push(root);
    
    while (!queue.empty()) {
        auto current = queue.front();
        queue.pop();
        
        if (current->getType() == HBNodeType::CONTOUR) {
            return current;
        }
        
        if (current->getLeftChild()) queue.push(current->getLeftChild());
        if (current->getRightChild()) queue.push(current->getRightChild());
    }
    
    return nullptr;
}

/**
 * Finds the leftmost skewed child of a node
 */
shared_ptr<HBStarTreeNode> HBStarTree::findLeftmostSkewedChild(shared_ptr<HBStarTreeNode> node) const {
    if (!node) return nullptr;
    
    auto current = node;
    while (current->getLeftChild()) {
        current = current->getLeftChild();
    }
    
    return current;
}

/**
 * Updates contour nodes after changing the ASF-B*-tree of a symmetry group
 */
void HBStarTree::updateContourNodes() {
    // Process each hierarchy node
    for (const auto& pair : symmetryGroupNodes) {
        auto hierarchyNode = pair.second;
        auto asfTree = hierarchyNode->getASFTree();
        
        if (!asfTree) continue;
        
        // Get the contours of the symmetry island
        auto contours = asfTree->getContours();
        auto horizontalContour = contours.first;
        
        // Get the horizontal contour segments
        auto segments = horizontalContour->getSegments();
        
        // Clear existing contour nodes
        vector<shared_ptr<HBStarTreeNode>> existingContourNodes;
        queue<shared_ptr<HBStarTreeNode>> queue;
        
        if (hierarchyNode->getRightChild()) {
            queue.push(hierarchyNode->getRightChild());
        }
        
        while (!queue.empty()) {
            auto current = queue.front();
            queue.pop();
            
            if (current->getType() == HBNodeType::CONTOUR) {
                existingContourNodes.push_back(current);
                
                if (current->getLeftChild()) {
                    queue.push(current->getLeftChild());
                }
                if (current->getRightChild()) {
                    queue.push(current->getRightChild());
                }
            }
        }
        
        // Create new contour nodes
        vector<shared_ptr<HBStarTreeNode>> newContourNodes;
        for (size_t i = 0; i < segments.size(); ++i) {
            auto contourNode = make_shared<HBStarTreeNode>(HBNodeType::CONTOUR, 
                                                             pair.first + "_contour_" + to_string(i));
            contourNode->setContour(segments[i].start, segments[i].height, segments[i].end, segments[i].height);
            newContourNodes.push_back(contourNode);
        }
        
        // Connect contour nodes
        if (!newContourNodes.empty()) {
            // Connect the first contour node to the hierarchy node
            hierarchyNode->setRightChild(newContourNodes[0]);
            newContourNodes[0]->setParent(hierarchyNode);
            
            // Connect the rest of the contour nodes
            for (size_t i = 0; i < newContourNodes.size() - 1; ++i) {
                newContourNodes[i]->setLeftChild(newContourNodes[i + 1]);
                newContourNodes[i + 1]->setParent(newContourNodes[i]);
            }
        }
        
        // Find dangling nodes - nodes whose parents were contour nodes that no longer exist
        vector<shared_ptr<HBStarTreeNode>> danglingNodes;
        for (const auto& oldContourNode : existingContourNodes) {
            if (oldContourNode->getRightChild()) {
                danglingNodes.push_back(oldContourNode->getRightChild());
            }
        }
        
        // Reassign dangling nodes
        for (const auto& danglingNode : danglingNodes) {
            // Find the nearest contour node
            auto nearestContourNode = findNearestContourNode(danglingNode);
            
            if (nearestContourNode) {
                if (!nearestContourNode->getRightChild()) {
                    // Attach directly as right child
                    nearestContourNode->setRightChild(danglingNode);
                    danglingNode->setParent(nearestContourNode);
                } else {
                    // Find the leftmost skewed child
                    auto leftmostSkewedChild = findLeftmostSkewedChild(nearestContourNode->getRightChild());
                    
                    leftmostSkewedChild->setLeftChild(danglingNode);
                    danglingNode->setParent(leftmostSkewedChild);
                }
            }
        }
    }
}

/**
 * Handles dangling nodes after tree modifications
 */
void HBStarTree::handleDanglingNodes() {
    // Similar to updateContourNodes but more general
    // Since this is complex and depends on implementation details, I'm leaving it blank for now
}

/**
 * Validates that all symmetry islands are placed correctly
 */
bool HBStarTree::validateSymmetryIslandPlacement() const {
    // Check each symmetry group
    for (const auto& group : symmetryGroups) {
        auto it = symmetryGroupNodes.find(group->getName());
        if (it == symmetryGroupNodes.end()) continue;
        
        auto hierarchyNode = it->second;
        auto asfTree = hierarchyNode->getASFTree();
        
        if (!asfTree || !asfTree->isSymmetricFeasible()) {
            return false;
        }
    }
    
    return true;
}

/**
 * Calculates the coordinates of all modules by packing the HB*-tree
 */
bool HBStarTree::pack() {
    if (!root) return false;
    
    // Reset contours
    horizontalContour->clear();
    verticalContour->clear();
    
    // Initialize horizontal contour with a segment at y=0
    horizontalContour->addSegment(0, numeric_limits<int>::max(), 0);
    
    // Initialize vertical contour with a segment at x=0
    verticalContour->addSegment(0, numeric_limits<int>::max(), 0);
    
    // Pack the tree in pre-order (DFS)
    stack<shared_ptr<HBStarTreeNode>> nodeStack;
    nodeStack.push(root);
    
    int maxX = 0, maxY = 0;
    
    while (!nodeStack.empty()) {
        auto currentNode = nodeStack.top();
        nodeStack.pop();
        
        if (!currentNode) continue;
        
        switch (currentNode->getType()) {
            case HBNodeType::MODULE: {
                // Pack a regular module
                const string& moduleName = currentNode->getModuleName();
                auto module = modules[moduleName];
                
                if (!module) continue;
                
                int x = 0, y = 0;
                
                // Calculate x-coordinate based on B*-tree rules
                if (currentNode->getParent()) {
                    if (currentNode->getParent()->getLeftChild() == currentNode) {
                        // Left child: place to the right of parent
                        if (currentNode->getParent()->getType() == HBNodeType::MODULE) {
                            auto parentModule = modules[currentNode->getParent()->getModuleName()];
                            if (parentModule) {
                                x = parentModule->getX() + parentModule->getWidth();
                            }
                        } else if (currentNode->getParent()->getType() == HBNodeType::HIERARCHY) {
                            // This is more complex - need to get the rightmost x-coordinate of the symmetry island
                            auto asfTree = currentNode->getParent()->getASFTree();
                            if (asfTree) {
                                // Use the symmetry axis position as a simple approximation
                                x = static_cast<int>(asfTree->getSymmetryAxisPosition());
                            }
                        } else if (currentNode->getParent()->getType() == HBNodeType::CONTOUR) {
                            // Place to the right of the contour segment
                            int x1, y1, x2, y2;
                            currentNode->getParent()->getContour(x1, y1, x2, y2);
                            x = x2;
                        }
                    } else {
                        // Right child: same x-coordinate as parent
                        if (currentNode->getParent()->getType() == HBNodeType::MODULE) {
                            auto parentModule = modules[currentNode->getParent()->getModuleName()];
                            if (parentModule) {
                                x = parentModule->getX();
                            }
                        } else if (currentNode->getParent()->getType() == HBNodeType::HIERARCHY) {
                            // Place at the left boundary of the symmetry island
                            // For simplicity, use x=0 for now
                            x = 0;
                        } else if (currentNode->getParent()->getType() == HBNodeType::CONTOUR) {
                            // Place at the left boundary of the contour segment
                            int x1, y1, x2, y2;
                            currentNode->getParent()->getContour(x1, y1, x2, y2);
                            x = x1;
                        }
                    }
                }
                
                // Calculate y-coordinate using the horizontal contour
                y = horizontalContour->getHeight(x, x + module->getWidth());
                
                // Set the module's position
                module->setPosition(x, y);
                
                // Update contours
                horizontalContour->addSegment(x, x + module->getWidth(), y + module->getHeight());
                verticalContour->addSegment(y, y + module->getHeight(), x + module->getWidth());
                
                // Update maximum coordinates
                maxX = max(maxX, x + module->getWidth());
                maxY = max(maxY, y + module->getHeight());
                
                break;
            }
            case HBNodeType::HIERARCHY: {
                // Pack a symmetry island
                auto asfTree = currentNode->getASFTree();
                if (!asfTree) continue;
                
                // Pack the ASF-B*-tree
                asfTree->pack();
                
                // Get the bounding rectangle of the symmetry island
                int minX = numeric_limits<int>::max();
                int minY = numeric_limits<int>::max();
                int symMaxX = 0;
                int symMaxY = 0;
                
                for (const auto& pair : asfTree->getModules()) {
                    const auto& module = pair.second;
                    
                    minX = min(minX, module->getX());
                    minY = min(minY, module->getY());
                    symMaxX = max(symMaxX, module->getX() + module->getWidth());
                    symMaxY = max(symMaxY, module->getY() + module->getHeight());
                }
                
                // Calculate the position for the symmetry island
                int x = 0, y = 0;
                
                // Calculate x-coordinate based on B*-tree rules
                if (currentNode->getParent()) {
                    if (currentNode->getParent()->getLeftChild() == currentNode) {
                        // Left child: place to the right of parent
                        if (currentNode->getParent()->getType() == HBNodeType::MODULE) {
                            auto parentModule = modules[currentNode->getParent()->getModuleName()];
                            if (parentModule) {
                                x = parentModule->getX() + parentModule->getWidth();
                            }
                        } else if (currentNode->getParent()->getType() == HBNodeType::HIERARCHY) {
                            // This is more complex - need to get the rightmost x-coordinate of the parent symmetry island
                            auto parentAsfTree = currentNode->getParent()->getASFTree();
                            if (parentAsfTree) {
                                // Use the symmetry axis position as a simple approximation
                                x = static_cast<int>(parentAsfTree->getSymmetryAxisPosition());
                            }
                        } else if (currentNode->getParent()->getType() == HBNodeType::CONTOUR) {
                            // Place to the right of the contour segment
                            int x1, y1, x2, y2;
                            currentNode->getParent()->getContour(x1, y1, x2, y2);
                            x = x2;
                        }
                    } else {
                        // Right child: same x-coordinate as parent
                        if (currentNode->getParent()->getType() == HBNodeType::MODULE) {
                            auto parentModule = modules[currentNode->getParent()->getModuleName()];
                            if (parentModule) {
                                x = parentModule->getX();
                            }
                        } else if (currentNode->getParent()->getType() == HBNodeType::HIERARCHY) {
                            // Place at the left boundary of the parent symmetry island
                            // For simplicity, use x=0 for now
                            x = 0;
                        } else if (currentNode->getParent()->getType() == HBNodeType::CONTOUR) {
                            // Place at the left boundary of the contour segment
                            int x1, y1, x2, y2;
                            currentNode->getParent()->getContour(x1, y1, x2, y2);
                            x = x1;
                        }
                    }
                }
                
                // Calculate y-coordinate using the horizontal contour
                y = horizontalContour->getHeight(x, x + (symMaxX - minX));
                
                // Shift all modules in the symmetry island
                int deltaX = x - minX;
                int deltaY = y - minY;
                
                for (const auto& pair : asfTree->getModules()) {
                    const auto& module = pair.second;
                    module->setPosition(module->getX() + deltaX, module->getY() + deltaY);
                }
                
                // Update contours - this is complex for rectilinear symmetry islands
                // For simplicity, use a rectangular approximation
                horizontalContour->addSegment(x, x + (symMaxX - minX), y + (symMaxY - minY));
                verticalContour->addSegment(y, y + (symMaxY - minY), x + (symMaxX - minX));
                
                // Update maximum coordinates
                maxX = max(maxX, x + (symMaxX - minX));
                maxY = max(maxY, y + (symMaxY - minY));
                
                break;
            }
            case HBNodeType::CONTOUR: {
                // Contour nodes don't need to be packed
                break;
            }
        }
        
        // Add children to the stack in reverse order (for pre-order traversal)
        if (currentNode->getRightChild()) {
            nodeStack.push(currentNode->getRightChild());
        }
        if (currentNode->getLeftChild()) {
            nodeStack.push(currentNode->getLeftChild());
        }
    }
    
    // Calculate total area
    totalArea = maxX * maxY;
    
    // Update contour nodes
    updateContourNodes();
    
    isPacked = true;
    
    return true;
}

/**
 * Performs a rotation operation on a module
 */
bool HBStarTree::rotateModule(const string& moduleName) {
    // Check if the module exists
    auto it = modules.find(moduleName);
    if (it == modules.end()) return false;
    
    auto module = it->second;
    
    // Check if the module is in a symmetry group
    bool inSymmetryGroup = false;
    shared_ptr<SymmetryGroup> group = nullptr;
    
    for (const auto& g : symmetryGroups) {
        // Check symmetry pairs
        for (const auto& pair : g->getSymmetryPairs()) {
            if (pair.first == moduleName || pair.second == moduleName) {
                inSymmetryGroup = true;
                group = g;
                break;
            }
        }
        
        // Check self-symmetric modules
        if (!inSymmetryGroup) {
            for (const auto& name : g->getSelfSymmetric()) {
                if (name == moduleName) {
                    inSymmetryGroup = true;
                    group = g;
                    break;
                }
            }
        }
        
        if (inSymmetryGroup) break;
    }
    
    // If the module is in a symmetry group, use the ASF-B*-tree to rotate it
    if (inSymmetryGroup && group) {
        auto it = symmetryGroupNodes.find(group->getName());
        if (it == symmetryGroupNodes.end()) return false;
        
        auto hierarchyNode = it->second;
        auto asfTree = hierarchyNode->getASFTree();
        
        if (!asfTree) return false;
        
        return asfTree->rotateModule(moduleName);
    }
    
    // Otherwise, just rotate the module directly
    module->rotate();
    
    // Since the module's dimensions have changed, the tree needs to be repacked
    if (isPacked) {
        pack();
    }
    
    return true;
}

/**
 * Moves a node to a new position in the tree
 */
bool HBStarTree::moveNode(const string& nodeName, 
                          const string& newParentName, 
                          bool asLeftChild) {
    // Find the nodes
    shared_ptr<HBStarTreeNode> node = nullptr;
    shared_ptr<HBStarTreeNode> newParent = nullptr;
    
    // Check if nodeName is a module name
    auto moduleIt = moduleNodes.find(nodeName);
    if (moduleIt != moduleNodes.end()) {
        node = moduleIt->second;
    }
    
    // Check if nodeName is a symmetry group name
    if (!node) {
        auto symGroupIt = symmetryGroupNodes.find(nodeName);
        if (symGroupIt != symmetryGroupNodes.end()) {
            node = symGroupIt->second;
        }
    }
    
    // Check if newParentName is a module name
    auto parentModuleIt = moduleNodes.find(newParentName);
    if (parentModuleIt != moduleNodes.end()) {
        newParent = parentModuleIt->second;
    }
    
    // Check if newParentName is a symmetry group name
    if (!newParent) {
        auto parentSymGroupIt = symmetryGroupNodes.find(newParentName);
        if (parentSymGroupIt != symmetryGroupNodes.end()) {
            newParent = parentSymGroupIt->second;
        }
    }
    
    if (!node || !newParent) return false;
    
    // Remove the node from its current parent
    auto oldParent = node->getParent();
    if (oldParent) {
        if (oldParent->getLeftChild() == node) {
            oldParent->setLeftChild(nullptr);
        } else if (oldParent->getRightChild() == node) {
            oldParent->setRightChild(nullptr);
        }
    } else if (node == root) {
        // The node is the root - find a new root
        if (node->getLeftChild()) {
            root = node->getLeftChild();
        } else if (node->getRightChild()) {
            root = node->getRightChild();
        } else {
            // No children - the tree will be empty after removal
            root = nullptr;
        }
    }
    
    // Add the node to its new parent
    node->setParent(newParent);
    if (asLeftChild) {
        // If there's already a left child, handle it
        auto existingChild = newParent->getLeftChild();
        if (existingChild) {
            // Try to find a place for the existing child
            if (!node->getLeftChild()) {
                node->setLeftChild(existingChild);
                existingChild->setParent(node);
            } else if (!node->getRightChild()) {
                node->setRightChild(existingChild);
                existingChild->setParent(node);
            } else {
                // Both children slots are taken - find another place
                // This is a simplification - in practice, more sophisticated handling is needed
                auto current = node->getLeftChild();
                while (current->getLeftChild()) {
                    current = current->getLeftChild();
                }
                current->setLeftChild(existingChild);
                existingChild->setParent(current);
            }
        }
        newParent->setLeftChild(node);
    } else {
        // If there's already a right child, handle it
        auto existingChild = newParent->getRightChild();
        if (existingChild) {
            // Try to find a place for the existing child
            if (!node->getLeftChild()) {
                node->setLeftChild(existingChild);
                existingChild->setParent(node);
            } else if (!node->getRightChild()) {
                node->setRightChild(existingChild);
                existingChild->setParent(node);
            } else {
                // Both children slots are taken - find another place
                auto current = node->getRightChild();
                while (current->getRightChild()) {
                    current = current->getRightChild();
                }
                current->setRightChild(existingChild);
                existingChild->setParent(current);
            }
        }
        newParent->setRightChild(node);
    }
    
    // Since the tree structure has changed, it needs to be repacked
    if (isPacked) {
        pack();
    }
    
    return true;
}

/**
 * Swaps two nodes in the tree
 */
bool HBStarTree::swapNodes(const string& nodeName1, const string& nodeName2) {
    // Find the nodes
    shared_ptr<HBStarTreeNode> node1 = nullptr;
    shared_ptr<HBStarTreeNode> node2 = nullptr;
    
    // Check if nodeName1 is a module name
    auto moduleIt1 = moduleNodes.find(nodeName1);
    if (moduleIt1 != moduleNodes.end()) {
        node1 = moduleIt1->second;
    }
    
    // Check if nodeName1 is a symmetry group name
    if (!node1) {
        auto symGroupIt1 = symmetryGroupNodes.find(nodeName1);
        if (symGroupIt1 != symmetryGroupNodes.end()) {
            node1 = symGroupIt1->second;
        }
    }
    
    // Check if nodeName2 is a module name
    auto moduleIt2 = moduleNodes.find(nodeName2);
    if (moduleIt2 != moduleNodes.end()) {
        node2 = moduleIt2->second;
    }
    
    // Check if nodeName2 is a symmetry group name
    if (!node2) {
        auto symGroupIt2 = symmetryGroupNodes.find(nodeName2);
        if (symGroupIt2 != symmetryGroupNodes.end()) {
            node2 = symGroupIt2->second;
        }
    }
    
    if (!node1 || !node2) return false;
    
    // Get parents and positions
    auto parent1 = node1->getParent();
    auto parent2 = node2->getParent();
    
    bool isLeftChild1 = parent1 && parent1->getLeftChild() == node1;
    bool isLeftChild2 = parent2 && parent2->getLeftChild() == node2;
    
    // Special case: node2 is a child of node1
    if (node1->getLeftChild() == node2 || node1->getRightChild() == node2) {
        // Detach node2 from node1
        if (node1->getLeftChild() == node2) {
            node1->setLeftChild(nullptr);
        } else {
            node1->setRightChild(nullptr);
        }
        
        // Detach node1 from its parent
        if (parent1) {
            if (isLeftChild1) {
                parent1->setLeftChild(nullptr);
            } else {
                parent1->setRightChild(nullptr);
            }
        }
        
        // Attach node1 as child of node2 in the same position
        if (node1->getLeftChild() == node2) {
            node2->setLeftChild(node1);
        } else {
            node2->setRightChild(node1);
        }
        node1->setParent(node2);
        
        // Attach node2 to node1's old parent
        if (parent1) {
            if (isLeftChild1) {
                parent1->setLeftChild(node2);
            } else {
                parent1->setRightChild(node2);
            }
            node2->setParent(parent1);
        } else {
            // node1 was the root
            root = node2;
            node2->setParent(nullptr);
        }
    }
    // Special case: node1 is a child of node2
    else if (node2->getLeftChild() == node1 || node2->getRightChild() == node1) {
        // Detach node1 from node2
        if (node2->getLeftChild() == node1) {
            node2->setLeftChild(nullptr);
        } else {
            node2->setRightChild(nullptr);
        }
        
        // Detach node2 from its parent
        if (parent2) {
            if (isLeftChild2) {
                parent2->setLeftChild(nullptr);
            } else {
                parent2->setRightChild(nullptr);
            }
        }
        
        // Attach node2 as child of node1 in the same position
        if (node2->getLeftChild() == node1) {
            node1->setLeftChild(node2);
        } else {
            node1->setRightChild(node2);
        }
        node2->setParent(node1);
        
        // Attach node1 to node2's old parent
        if (parent2) {
            if (isLeftChild2) {
                parent2->setLeftChild(node1);
            } else {
                parent2->setRightChild(node1);
            }
            node1->setParent(parent2);
        } else {
            // node2 was the root
            root = node1;
            node1->setParent(nullptr);
        }
    }
    // General case: nodes are not directly related
    else {
        // Detach nodes from parents
        if (parent1) {
            if (isLeftChild1) {
                parent1->setLeftChild(nullptr);
            } else {
                parent1->setRightChild(nullptr);
            }
        }
        
        if (parent2) {
            if (isLeftChild2) {
                parent2->setLeftChild(nullptr);
            } else {
                parent2->setRightChild(nullptr);
            }
        }
        
        // Swap children
        auto leftChild1 = node1->getLeftChild();
        auto rightChild1 = node1->getRightChild();
        auto leftChild2 = node2->getLeftChild();
        auto rightChild2 = node2->getRightChild();
        
        // Set children for node1
        node1->setLeftChild(leftChild2);
        node1->setRightChild(rightChild2);
        if (leftChild2) leftChild2->setParent(node1);
        if (rightChild2) rightChild2->setParent(node1);
        
        // Set children for node2
        node2->setLeftChild(leftChild1);
        node2->setRightChild(rightChild1);
        if (leftChild1) leftChild1->setParent(node2);
        if (rightChild1) rightChild1->setParent(node2);
        
        // Reattach nodes to opposite parents
        if (parent1) {
            if (isLeftChild1) {
                parent1->setLeftChild(node2);
            } else {
                parent1->setRightChild(node2);
            }
            node2->setParent(parent1);
        } else {
            // node1 was the root
            root = node2;
            node2->setParent(nullptr);
        }
        
        if (parent2) {
            if (isLeftChild2) {
                parent2->setLeftChild(node1);
            } else {
                parent2->setRightChild(node1);
            }
            node1->setParent(parent2);
        } else {
            // node2 was the root
            root = node1;
            node1->setParent(nullptr);
        }
    }
    
    // Since the tree structure has changed, it needs to be repacked
    if (isPacked) {
        pack();
    }
    
    return true;
}

/**
 * Changes the representative of a symmetry pair in a symmetry group
 */
bool HBStarTree::changeRepresentative(const string& symmetryGroupName, 
                                     const string& moduleName) {
    // Find the symmetry group
    auto it = symmetryGroupNodes.find(symmetryGroupName);
    if (it == symmetryGroupNodes.end()) return false;
    
    auto hierarchyNode = it->second;
    auto asfTree = hierarchyNode->getASFTree();
    
    if (!asfTree) return false;
    
    // Change the representative
    bool success = asfTree->changeRepresentative(moduleName);
    
    // Since the symmetry island might have changed, the tree needs to be repacked
    if (success && isPacked) {
        pack();
    }
    
    return success;
}

/**
 * Converts the symmetry type of a symmetry group
 */
bool HBStarTree::convertSymmetryType(const string& symmetryGroupName) {
    // Find the symmetry group
    auto it = symmetryGroupNodes.find(symmetryGroupName);
    if (it == symmetryGroupNodes.end()) return false;
    
    auto hierarchyNode = it->second;
    auto asfTree = hierarchyNode->getASFTree();
    
    if (!asfTree) return false;
    
    // Convert the symmetry type
    bool success = asfTree->convertSymmetryType();
    
    // Since the symmetry island might have changed, the tree needs to be repacked
    if (success && isPacked) {
        pack();
    }
    
    return success;
}

/**
 * Returns the total area of the placement
 */
int HBStarTree::getArea() const {
    return totalArea;
}

/**
 * Returns the total wire length of the placement
 */
int HBStarTree::getWireLength() const {
    // Wire length calculation depends on netlist information, return 0 for now
    return 0;
}

/**
 * Gets the root node of the HB*-tree
 */
shared_ptr<HBStarTreeNode> HBStarTree::getRoot() const {
    return root;
}

/**
 * Gets all modules in the design
 */
const map<string, shared_ptr<Module>>& HBStarTree::getModules() const {
    return modules;
}

/**
 * Gets all symmetry groups in the design
 */
const vector<shared_ptr<SymmetryGroup>>& HBStarTree::getSymmetryGroups() const {
    return symmetryGroups;
}

/**
 * Gets the module node with the given name
 */
shared_ptr<HBStarTreeNode> HBStarTree::getModuleNode(const string& moduleName) const {
    auto it = moduleNodes.find(moduleName);
    if (it == moduleNodes.end()) return nullptr;
    
    return it->second;
}

/**
 * Gets the symmetry group node with the given name
 */
shared_ptr<HBStarTreeNode> HBStarTree::getSymmetryGroupNode(const string& symmetryGroupName) const {
    auto it = symmetryGroupNodes.find(symmetryGroupName);
    if (it == symmetryGroupNodes.end()) return nullptr;
    
    return it->second;
}

/**
 * Creates a deep copy of this HB*-tree
 */
shared_ptr<HBStarTree> HBStarTree::clone() const {
    auto clone = make_shared<HBStarTree>();
    
    // Copy modules
    for (const auto& pair : modules) {
        auto moduleCopy = make_shared<Module>(*pair.second);
        clone->modules[pair.first] = moduleCopy;
    }
    
    // Copy symmetry groups
    for (const auto& group : symmetryGroups) {
        // Deep copy of symmetry group
        auto groupCopy = make_shared<SymmetryGroup>(*group);
        clone->symmetryGroups.push_back(groupCopy);
    }
    
    // Reconstruct initial tree
    clone->constructInitialTree();
    
    // Copy packed state
    clone->isPacked = isPacked;
    clone->totalArea = totalArea;
    
    // Copy contours
    clone->horizontalContour = make_shared<Contour>(*horizontalContour);
    clone->verticalContour = make_shared<Contour>(*verticalContour);
    
    return clone;
}