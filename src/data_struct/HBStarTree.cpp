/**
 * HBStarTree.cpp
 * 
 * Implementation of the HBStarTree class for analog placement with symmetry constraints.
 * The HB*-tree is a hierarchical framework that can simultaneously optimize the placement
 * with both symmetry islands and non-symmetric modules.
 * 
 * This version includes incremental packing for improved performance.
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
      packed(false) {
}

HBStarTree::~HBStarTree() {
    clearTree();
}

/* Add a module to the tree */
void HBStarTree::addModule(shared_ptr<Module> module) {
    if (!module) return;
    
    // Add the module to module map
    modules[module->getName()] = module;
    
    // Mark as needing repacking
    packed = false;
}

/* Add a symmetry group to the tree */
void HBStarTree::addSymmetryGroup(shared_ptr<SymmetryGroup> group) {
    if (!group) return;
    
    // Add the symmetry group to list
    symmetryGroups.push_back(group);
    
    // Mark as needing repacking
    packed = false;
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
    packed = false;
    clearAffectedCache();
}

/* Constructs an initial HB*-tree */
void HBStarTree::constructInitialTree() {
    // Clear any existing tree
    clearTree();
    
    // First, construct symmetry islands for each symmetry group
    constructSymmetryIslands();
    
    // Then, construct the initial tree structure
    constructInitialTreeStructure();
    
    // Mark as needing full packing
    packed = false;
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
        
        // Skip if this symmetry group is not affected
        if (affectedModules.find(pair.first) == affectedModules.end() && 
            affectedNodes.find(hierarchyNode) == affectedNodes.end()) {
            continue;
        }
        
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
 * Pack a single node in the HB*-tree
 */
bool HBStarTree::packNode(shared_ptr<HBStarTreeNode> node) {
    if (!node) return false;
    
    switch (node->getType()) {
        case HBNodeType::MODULE: {
            const string& moduleName = node->getModuleName();
            auto module = modules[moduleName];
            
            if (!module) return false;
            
            int x = 0, y = 0;
            
            // Calculate x-coordinate based on B*-tree rules
            if (node->getParent()) {
                if (node->getParent()->getLeftChild() == node) {
                    // Left child: place to the right of parent
                    if (node->getParent()->getType() == HBNodeType::MODULE) {
                        auto parentModule = modules[node->getParent()->getModuleName()];
                        if (parentModule) {
                            x = parentModule->getX() + parentModule->getWidth();
                        }
                    } else if (node->getParent()->getType() == HBNodeType::HIERARCHY) {
                        // This is more complex - need to get the rightmost x-coordinate of the symmetry island
                        auto asfTree = node->getParent()->getASFTree();
                        if (asfTree) {
                            // Use the symmetry axis position as a simple approximation
                            x = static_cast<int>(asfTree->getSymmetryAxisPosition());
                        }
                    } else if (node->getParent()->getType() == HBNodeType::CONTOUR) {
                        // Place to the right of the contour segment
                        int x1, y1, x2, y2;
                        node->getParent()->getContour(x1, y1, x2, y2);
                        x = x2;
                    }
                } else {
                    // Right child: same x-coordinate as parent
                    if (node->getParent()->getType() == HBNodeType::MODULE) {
                        auto parentModule = modules[node->getParent()->getModuleName()];
                        if (parentModule) {
                            x = parentModule->getX();
                        }
                    } else if (node->getParent()->getType() == HBNodeType::HIERARCHY) {
                        // Place at the left boundary of the symmetry island
                        // For simplicity, use x=0 for now
                        x = 0;
                    } else if (node->getParent()->getType() == HBNodeType::CONTOUR) {
                        // Place at the left boundary of the contour segment
                        int x1, y1, x2, y2;
                        node->getParent()->getContour(x1, y1, x2, y2);
                        x = x1;
                    }
                }
            }
            
            // Calculate y-coordinate using the horizontal contour
            y = horizontalContour->getHeight(x, x + module->getWidth());
            
            // Set the module's position
            module->setPosition(x, y);
            
            // Update contours
            updateContourForModule(module);
            
            return true;
        }
        case HBNodeType::HIERARCHY: {
            return packSymmetryIsland(node);
        }
        case HBNodeType::CONTOUR: {
            // Contour nodes don't need to be packed themselves
            return true;
        }
        default:
            return false;
    }
}

/**
 * Updates contour with a module
 */
bool HBStarTree::updateContourForModule(const shared_ptr<Module>& module) {
    if (!module) return false;
    
    int x = module->getX();
    int y = module->getY();
    int width = module->getWidth();
    int height = module->getHeight();
    
    // Update horizontal contour
    horizontalContour->addSegment(x, x + width, y + height);
    
    // Update vertical contour
    verticalContour->addSegment(y, y + height, x + width);
    
    return true;
}

/**
 * Pack a symmetry island
 */
bool HBStarTree::packSymmetryIsland(shared_ptr<HBStarTreeNode> hierarchyNode) {
    if (!hierarchyNode || hierarchyNode->getType() != HBNodeType::HIERARCHY) return false;
    
    auto asfTree = hierarchyNode->getASFTree();
    if (!asfTree) return false;
    
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
    if (hierarchyNode->getParent()) {
        if (hierarchyNode->getParent()->getLeftChild() == hierarchyNode) {
            // Left child: place to the right of parent
            if (hierarchyNode->getParent()->getType() == HBNodeType::MODULE) {
                auto parentModule = modules[hierarchyNode->getParent()->getModuleName()];
                if (parentModule) {
                    x = parentModule->getX() + parentModule->getWidth();
                }
            } else if (hierarchyNode->getParent()->getType() == HBNodeType::HIERARCHY) {
                // This is more complex - need to get the rightmost x-coordinate of the parent symmetry island
                auto parentAsfTree = hierarchyNode->getParent()->getASFTree();
                if (parentAsfTree) {
                    // Use the symmetry axis position as a simple approximation
                    x = static_cast<int>(parentAsfTree->getSymmetryAxisPosition());
                }
            } else if (hierarchyNode->getParent()->getType() == HBNodeType::CONTOUR) {
                // Place to the right of the contour segment
                int x1, y1, x2, y2;
                hierarchyNode->getParent()->getContour(x1, y1, x2, y2);
                x = x2;
            }
        } else {
            // Right child: same x-coordinate as parent
            if (hierarchyNode->getParent()->getType() == HBNodeType::MODULE) {
                auto parentModule = modules[hierarchyNode->getParent()->getModuleName()];
                if (parentModule) {
                    x = parentModule->getX();
                }
            } else if (hierarchyNode->getParent()->getType() == HBNodeType::HIERARCHY) {
                // Place at the left boundary of the parent symmetry island
                // For simplicity, use x=0 for now
                x = 0;
            } else if (hierarchyNode->getParent()->getType() == HBNodeType::CONTOUR) {
                // Place at the left boundary of the contour segment
                int x1, y1, x2, y2;
                hierarchyNode->getParent()->getContour(x1, y1, x2, y2);
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
    
    return true;
}

/**
 * Pack a subtree of the HB*-tree
 */
bool HBStarTree::packSubtree(shared_ptr<HBStarTreeNode> node) {
    if (!node) return false;
    
    // First pack the node itself
    if (!packNode(node)) return false;
    
    // Then pack its children
    if (node->getLeftChild()) {
        if (!packSubtree(node->getLeftChild())) return false;
    }
    if (node->getRightChild()) {
        if (!packSubtree(node->getRightChild())) return false;
    }
    
    return true;
}

/**
 * Mark a module as affected, needing repacking
 */
void HBStarTree::markModuleAffected(const string& moduleName) {
    affectedModules.insert(moduleName);
    
    // Mark the associated node as affected
    auto it = moduleNodes.find(moduleName);
    if (it != moduleNodes.end()) {
        affectedNodes.insert(it->second);
    }
    
    // Also mark any symmetry group containing this module
    for (const auto& group : symmetryGroups) {
        bool inGroup = false;
        
        // Check symmetry pairs
        for (const auto& pair : group->getSymmetryPairs()) {
            if (pair.first == moduleName || pair.second == moduleName) {
                inGroup = true;
                break;
            }
        }
        
        // Check self-symmetric modules
        if (!inGroup) {
            for (const auto& name : group->getSelfSymmetric()) {
                if (name == moduleName) {
                    inGroup = true;
                    break;
                }
            }
        }
        
        if (inGroup) {
            markSymmetryGroupAffected(group->getName());
        }
    }
    
    // Mark as needing repacking
    packed = false;
}

/**
 * Mark a symmetry group as affected, needing repacking
 */
void HBStarTree::markSymmetryGroupAffected(const string& symmetryGroupName) {
    affectedModules.insert(symmetryGroupName);
    
    // Mark the associated node as affected
    auto it = symmetryGroupNodes.find(symmetryGroupName);
    if (it != symmetryGroupNodes.end()) {
        affectedNodes.insert(it->second);
    }
    
    // Mark as needing repacking
    packed = false;
}

/**
 * Mark a node as affected, needing repacking
 */
void HBStarTree::markNodeAffected(shared_ptr<HBStarTreeNode> node) {
    if (!node) return;
    
    affectedNodes.insert(node);
    
    // If it's a module node, mark the module as affected
    if (node->getType() == HBNodeType::MODULE) {
        affectedModules.insert(node->getModuleName());
    }
    // If it's a hierarchy node, mark the symmetry group as affected
    else if (node->getType() == HBNodeType::HIERARCHY) {
        affectedModules.insert(node->getName());
    }
    
    // Mark as needing repacking
    packed = false;
}

/**
 * Mark an affected subtree
 */
void HBStarTree::markAffectedSubtree(shared_ptr<HBStarTreeNode> node) {
    if (!node) return;
    
    markNodeAffected(node);
    
    // Mark all descendants
    if (node->getLeftChild()) markAffectedSubtree(node->getLeftChild());
    if (node->getRightChild()) markAffectedSubtree(node->getRightChild());
}

/**
 * Mark affected modules for a module
 */
void HBStarTree::markAffectedModules(const string& moduleName) {
    // Mark the module itself
    markModuleAffected(moduleName);
    
    // Mark all modules that are affected by this module's position change
    auto it = moduleNodes.find(moduleName);
    if (it != moduleNodes.end()) {
        auto node = it->second;
        
        // Mark descendants
        markAffectedSubtree(node);
        
        // Mark siblings to the right
        auto parent = node->getParent();
        if (parent) {
            auto sibling = parent->getRightChild();
            if (sibling && sibling != node) {
                markAffectedSubtree(sibling);
            }
        }
    }
}

/**
 * Clear the affected cache
 */
void HBStarTree::clearAffectedCache() {
    affectedModules.clear();
    affectedNodes.clear();
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
    
    // Clear affected cache
    clearAffectedCache();
    
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
                packNode(currentNode);
                
                // Update maximum coordinates
                const string& moduleName = currentNode->getModuleName();
                auto module = modules[moduleName];
                if (module) {
                    maxX = max(maxX, module->getX() + module->getWidth());
                    maxY = max(maxY, module->getY() + module->getHeight());
                }
                break;
            }
            case HBNodeType::HIERARCHY: {
                // Pack a symmetry island
                packSymmetryIsland(currentNode);
                
                // Update maximum coordinates using the coordinates of the modules in the symmetry island
                auto asfTree = currentNode->getASFTree();
                if (asfTree) {
                    for (const auto& pair : asfTree->getModules()) {
                        const auto& module = pair.second;
                        maxX = max(maxX, module->getX() + module->getWidth());
                        maxY = max(maxY, module->getY() + module->getHeight());
                    }
                }
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
    
    packed = true;
    
    return true;
}

/**
 * Incrementally updates the packing after a perturbation
 */
bool HBStarTree::incrementalPack() {
    if (!root) return false;
    
    // If no affected modules or nodes, return immediately
    if (affectedModules.empty() && affectedNodes.empty()) {
        if (!packed) {
            // If not packed, do a full pack
            return pack();
        }
        return true;
    }
    
    // Get the maximum coordinates from the current placement
    int maxX = 0, maxY = 0;
    
    // Find root nodes of affected subtrees
    unordered_set<shared_ptr<HBStarTreeNode>> rootsToRepack;
    
    // First, identify all the affected nodes that need repacking
    for (auto node : affectedNodes) {
        // Find the highest ancestor that needs repacking
        shared_ptr<HBStarTreeNode> current = node;
        shared_ptr<HBStarTreeNode> highest = node;
        
        while (current->getParent()) {
            current = current->getParent();
            
            // If parent is affected or current is right child, update highest
            if (affectedNodes.find(current) != affectedNodes.end() || 
                current->getRightChild() == highest) {
                highest = current;
            }
        }
        
        rootsToRepack.insert(highest);
    }
    
    // Repack each affected subtree
    for (auto node : rootsToRepack) {
        packSubtree(node);
    }
    
    // Update maximum coordinates
    for (const auto& pair : modules) {
        const auto& module = pair.second;
        maxX = max(maxX, module->getX() + module->getWidth());
        maxY = max(maxY, module->getY() + module->getHeight());
    }
    
    // Calculate total area
    totalArea = maxX * maxY;
    
    // Update contour nodes
    updateContourNodes();
    
    // Clear affected cache
    clearAffectedCache();
    
    packed = true;
    
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
        
        // Mark the symmetry group as affected
        markSymmetryGroupAffected(group->getName());
        
        return asfTree->rotateModule(moduleName);
    }
    
    // Otherwise, just rotate the module directly
    module->rotate();
    
    // Mark the module as affected
    markModuleAffected(moduleName);
    
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
    
    // Mark both nodes as affected
    markNodeAffected(node);
    markNodeAffected(newParent);
    
    // Remove the node from its current parent
    auto oldParent = node->getParent();
    if (oldParent) {
        if (oldParent->getLeftChild() == node) {
            oldParent->setLeftChild(nullptr);
        } else if (oldParent->getRightChild() == node) {
            oldParent->setRightChild(nullptr);
        }
        
        // Mark the old parent as affected
        markNodeAffected(oldParent);
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
            // Mark the existing child as affected
            markNodeAffected(existingChild);
            
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
            // Mark the existing child as affected
            markNodeAffected(existingChild);
            
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
    
    // Mark both nodes as affected
    markNodeAffected(node1);
    markNodeAffected(node2);
    
    // Get parents and positions
    auto parent1 = node1->getParent();
    auto parent2 = node2->getParent();
    
    bool isLeftChild1 = parent1 && parent1->getLeftChild() == node1;
    bool isLeftChild2 = parent2 && parent2->getLeftChild() == node2;
    
    // Mark parents as affected
    if (parent1) markNodeAffected(parent1);
    if (parent2) markNodeAffected(parent2);
    
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
        
        // Mark children as affected
        if (leftChild1) markNodeAffected(leftChild1);
        if (rightChild1) markNodeAffected(rightChild1);
        if (leftChild2) markNodeAffected(leftChild2);
        if (rightChild2) markNodeAffected(rightChild2);
        
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
    
    // Mark the symmetry group as affected
    markSymmetryGroupAffected(symmetryGroupName);
    
    // Change the representative
    return asfTree->changeRepresentative(moduleName);
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
    
    // Mark the symmetry group as affected
    markSymmetryGroupAffected(symmetryGroupName);
    
    // Convert the symmetry type
    return asfTree->convertSymmetryType();
}

/**
 * Returns whether the tree has been packed
 */
bool HBStarTree::isPacked() const {
    return packed;
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
    clone->packed = packed;
    clone->totalArea = totalArea;
    
    // Copy contours
    clone->horizontalContour = make_shared<Contour>(*horizontalContour);
    clone->verticalContour = make_shared<Contour>(*verticalContour);
    
    return clone;
}