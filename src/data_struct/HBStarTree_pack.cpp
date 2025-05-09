/**
 * HBStarTree_Packing.cpp
 * 
 * Implementation of the packing and perturbation operations of the HBStarTree class
 * for analog placement with symmetry constraints.
 */

#include "HBStarTree.hpp"
#include <algorithm>
#include <queue>
#include <iostream>
#include <limits>
#include <stack>
using namespace std;

/**
 * Calculates the coordinates of all modules by packing the HB*-tree
 */
bool HBStarTree::pack() {
    if (!root) {
        std::cerr << "Cannot pack HBStarTree: root is null" << std::endl;
        return false;
    }
    
    // Reset contours
    horizontalContour->clear();
    verticalContour->clear();
    
    // Initialize horizontal contour with a segment at y=0
    horizontalContour->addSegment(0, numeric_limits<int>::max(), 0);
    
    // Initialize vertical contour with a segment at x=0
    verticalContour->addSegment(0, numeric_limits<int>::max(), 0);
    
    // Clear affected cache
    clearAffectedCache();
    
    // Pack all symmetry groups first
    for (auto& pair : symmetryGroupNodes) {
        auto hierarchyNode = pair.second;
        auto asfTree = hierarchyNode->getASFTree();
        
        if (asfTree) {
            // Pack the ASF-B*-tree to get initial module positions within the symmetry group
            if (!asfTree->pack()) {
                std::cerr << "Failed to pack symmetry group: " << pair.first << std::endl;
                return false; // Failed to pack a symmetry group
            }
        }
    }
    
    // IMPORTANT CHANGE: Use recursive packing for proper tree traversal
    bool packSuccess = packSubtree(root);
    
    if (!packSuccess) {
        std::cerr << "Failed to pack the HB*-tree" << std::endl;
        return false;
    }
    
    // Calculate maximum coordinates
    int maxX = 0, maxY = 0;
    
    for (const auto& pair : modules) {
        const auto& module = pair.second;
        maxX = max(maxX, module->getX() + module->getWidth());
        maxY = max(maxY, module->getY() + module->getHeight());
    }
    
    // For symmetry groups, check all modules in their ASF-B*-trees
    for (const auto& pair : symmetryGroupNodes) {
        auto asfTree = pair.second->getASFTree();
        if (asfTree) {
            for (const auto& modulePair : asfTree->getModules()) {
                const auto& module = modulePair.second;
                maxX = max(maxX, module->getX() + module->getWidth());
                maxY = max(maxY, module->getY() + module->getHeight());
            }
        }
    }
    
    // Calculate total area
    totalArea = maxX * maxY;
    
    // Update contour nodes
    updateContourNodes();
    
    // Validate symmetry constraints for each symmetry group
    bool allValid = true;
    for (const auto& symGroupPair : symmetryGroupNodes) {
        auto asfTree = symGroupPair.second->getASFTree();
        if (asfTree) {
            // Validate and fix symmetry constraints if needed
            asfTree->enforceSymmetryConstraints();
            
            // Check if the symmetry constraints are still valid after fixing
            if (!asfTree->checkSymmetryConstraints()) {
                std::cerr << "Symmetry constraints still violated for group: " 
                          << symGroupPair.first << " after fixing" << std::endl;
                allValid = false;
            }
            
            // Check for symmetry island validity
            if (!asfTree->isSymmetryIslandValid()) {
                std::cerr << "Invalid symmetry island for group: " 
                          << symGroupPair.first << std::endl;
                allValid = false;
            }
        }
    }
    
    // Check for overlaps after the symmetry enforcement
    if (hasOverlap()) {
        std::cerr << "Overlap detected in placement after packing" << std::endl;
        allValid = false;
    }
    
    if (!allValid) {
        // Attempt emergency recovery if constraints violated
        for (const auto& symGroupPair : symmetryGroupNodes) {
            auto asfTree = symGroupPair.second->getASFTree();
            if (asfTree && (!asfTree->checkSymmetryConstraints() || !asfTree->isSymmetryIslandValid())) {
                std::cout << "Applying emergency recovery for group: " << symGroupPair.first << std::endl;
                asfTree->emergencyRecovery();
            }
        }
        
        // Re-check after emergency recovery
        bool recoverySuccessful = true;
        for (const auto& symGroupPair : symmetryGroupNodes) {
            auto asfTree = symGroupPair.second->getASFTree();
            if (asfTree && (!asfTree->checkSymmetryConstraints() || !asfTree->isSymmetryIslandValid())) {
                recoverySuccessful = false;
                break;
            }
        }
        
        if (!recoverySuccessful) {
            std::cerr << "Emergency recovery failed to fix all issues" << std::endl;
            return false;
        }
    }
    
    packed = true;
    
    return true;
}

/**
 * Checks if the placement has any module overlaps
 */
bool HBStarTree::hasOverlap() const {
    // Check all pairs of modules for overlap
    std::vector<std::shared_ptr<Module>> allModules;
    
    // First, collect all modules (both regular and from symmetry groups)
    for (const auto& pair : modules) {
        allModules.push_back(pair.second);
    }
    
    // Also add modules from symmetry groups
    for (const auto& pair : symmetryGroupNodes) {
        const auto& node = pair.second;
        if (node && node->getType() == HBNodeType::HIERARCHY) {
            auto asfTree = node->getASFTree();
            if (asfTree) {
                for (const auto& modPair : asfTree->getModules()) {
                    allModules.push_back(modPair.second);
                }
            }
        }
    }
    
    // Check all pairs of modules for overlap
    for (size_t i = 0; i < allModules.size(); ++i) {
        const auto& module1 = allModules[i];
        
        for (size_t j = i + 1; j < allModules.size(); ++j) {
            const auto& module2 = allModules[j];
            
            // Skip if comparing a module with itself or comparing modules with the same name
            if (module1 == module2 || module1->getName() == module2->getName()) {
                continue;
            }
            
            if (module1->overlaps(*module2)) {
                std::cerr << "Overlap detected between modules: " 
                          << module1->getName() << " (" << module1->getX() << "," << module1->getY() 
                          << "," << module1->getWidth() << "," << module1->getHeight() << ") and " 
                          << module2->getName() << " (" << module2->getX() << "," << module2->getY()
                          << "," << module2->getWidth() << "," << module2->getHeight() << ")" << std::endl;
                return true;
            }
        }
    }
    
    return false;
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
 * Pack a single node in the HB*-tree
 */
bool HBStarTree::packNode(shared_ptr<HBStarTreeNode> node) {
    if (!node) return false;
    
    switch (node->getType()) {
        case HBNodeType::MODULE: {
            const string& moduleName = node->getModuleName();
            auto moduleIt = modules.find(moduleName);
            if (moduleIt == modules.end()) return false;
            
            auto module = moduleIt->second;
            
            int x = 0, y = 0;
            
            // Calculate x-coordinate based on B*-tree rules
            if (node->getParent()) {
                if (node->isLeftChild()) {
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
                            // Find the maximum x-coordinate of all modules in the symmetry island
                            int maxX = 0;
                            for (const auto& modulePair : asfTree->getModules()) {
                                const auto& symModule = modulePair.second;
                                maxX = max(maxX, symModule->getX() + symModule->getWidth());
                            }
                            x = maxX;
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
                        auto asfTree = node->getParent()->getASFTree();
                        if (asfTree) {
                            // Find the minimum x-coordinate of all modules in the symmetry island
                            int minX = numeric_limits<int>::max();
                            for (const auto& modulePair : asfTree->getModules()) {
                                const auto& symModule = modulePair.second;
                                minX = min(minX, symModule->getX());
                            }
                            x = minX;
                        }
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
            
            // Update horizontal contour
            horizontalContour->addSegment(x, x + module->getWidth(), y + module->getHeight());
            
            // Update vertical contour
            verticalContour->addSegment(y, y + module->getHeight(), x + module->getWidth());
            
            return true;
        }
        case HBNodeType::HIERARCHY:
        case HBNodeType::CONTOUR:
            // These nodes are handled separately
            return true;
        default:
            return false;
    }
}


/**
 * Pack a symmetry island
 */
/**
 * Pack a symmetry island
 */
bool HBStarTree::packSymmetryIsland(shared_ptr<HBStarTreeNode> hierarchyNode) {
    if (!hierarchyNode || hierarchyNode->getType() != HBNodeType::HIERARCHY) return false;
    
    auto asfTree = hierarchyNode->getASFTree();
    if (!asfTree) return false;
    
    // The ASF-B*-tree was already packed earlier, so now we need to position the entire island
    
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
        if (hierarchyNode->isLeftChild()) {
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
                    // Find the maximum x-coordinate of all modules in the parent symmetry island
                    int maxX = 0;
                    for (const auto& modulePair : parentAsfTree->getModules()) {
                        const auto& symModule = modulePair.second;
                        maxX = max(maxX, symModule->getX() + symModule->getWidth());
                    }
                    x = maxX;
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
                auto parentAsfTree = hierarchyNode->getParent()->getASFTree();
                if (parentAsfTree) {
                    // Find the minimum x-coordinate of all modules in the parent symmetry island
                    int minParentX = numeric_limits<int>::max();
                    for (const auto& modulePair : parentAsfTree->getModules()) {
                        const auto& symModule = modulePair.second;
                        minParentX = min(minParentX, symModule->getX());
                    }
                    x = minParentX;
                }
            } else if (hierarchyNode->getParent()->getType() == HBNodeType::CONTOUR) {
                // Place at the left boundary of the contour segment
                int x1, y1, x2, y2;
                hierarchyNode->getParent()->getContour(x1, y1, x2, y2);
                x = x1;
            }
        }
    }
    
    // Calculate y-coordinate using the horizontal contour
    int width = symMaxX - minX;
    y = horizontalContour->getHeight(x, x + width);
    
    // Check for y-coordinates from any below contour segments
    if (hierarchyNode->getParent() && hierarchyNode->getParent()->getType() == HBNodeType::CONTOUR) {
        int x1, y1, x2, y2;
        hierarchyNode->getParent()->getContour(x1, y1, x2, y2);
        y = max(y, y2); // Ensure placement above contour segment
    }
    
    // Calculate the offsets to shift the island
    int islandOffsetX = x - minX;
    int islandOffsetY = y - minY;
    
    // Update symmetry axis position if needed
    if (asfTree->getSymmetryGroup()->getType() == SymmetryType::VERTICAL) {
        // For vertical symmetry, update the axis position
        double oldAxis = asfTree->getSymmetryAxisPosition();
        asfTree->getSymmetryGroup()->setAxisPosition(oldAxis + islandOffsetX);
    } else { // HORIZONTAL
        // For horizontal symmetry, update the axis position
        double oldAxis = asfTree->getSymmetryAxisPosition();
        asfTree->getSymmetryGroup()->setAxisPosition(oldAxis + islandOffsetY);
    }
    
    // Shift all modules in the symmetry island
    for (const auto& pair : asfTree->getModules()) {
        auto module = pair.second;
        module->setPosition(module->getX() + islandOffsetX, module->getY() + islandOffsetY);
    }
    
    // IMPORTANT CHANGE: Do not add the bounding box as a simplification
    // Only add detailed contour segments from the symmetry island
    
    // Add the horizontal contour segments directly
    auto segments = asfTree->getContours().first->getSegments();
    for (const auto& segment : segments) {
        horizontalContour->addSegment(
            segment.start + islandOffsetX,
            segment.end + islandOffsetX,
            segment.height + islandOffsetY
        );
    }
    
    // Add the vertical contour segments
    const auto &vSegs = asfTree->getContours().second->getSegments();
    for (const auto &s : vSegs) {
        verticalContour->addSegment(
            s.start + islandOffsetY,
            s.end + islandOffsetY,
            s.height + islandOffsetX // x-coord lives in `height` field
        );
    }
    
    // After moving the island, re-validate symmetry
    asfTree->enforceSymmetryConstraints();
    
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
 * Incrementally updates the packing after a perturbation
 */
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
    
    // Repack each affected symmetry group first
    for (const auto& groupName : affectedModules) {
        auto symGroupIt = symmetryGroupNodes.find(groupName);
        if (symGroupIt != symmetryGroupNodes.end()) {
            auto asfTree = symGroupIt->second->getASFTree();
            if (asfTree) {
                if (!asfTree->pack()) {
                    return false; // Failed to pack the symmetry group
                }
                
                // Re-enforce symmetry constraints
                asfTree->enforceSymmetryConstraints();
            }
        }
    }
    
    // Repack each affected subtree
    for (auto node : rootsToRepack) {
        if (!packSubtree(node)) {
            return false; // Failed to pack subtree
        }
    }
    
    // Update maximum coordinates
    for (const auto& pair : modules) {
        const auto& module = pair.second;
        maxX = max(maxX, module->getX() + module->getWidth());
        maxY = max(maxY, module->getY() + module->getHeight());
    }
    
    // For symmetry groups, check all modules in their ASF-B*-trees
    for (const auto& pair : symmetryGroupNodes) {
        auto asfTree = pair.second->getASFTree();
        if (asfTree) {
            for (const auto& modulePair : asfTree->getModules()) {
                const auto& module = modulePair.second;
                maxX = max(maxX, module->getX() + module->getWidth());
                maxY = max(maxY, module->getY() + module->getHeight());
            }
        }
    }
    
    // Calculate total area
    totalArea = maxX * maxY;
    
    // Update contour nodes
    updateContourNodes();
    
    // Re-validate symmetry constraints for all affected symmetry groups
    for (const auto& groupName : affectedModules) {
        auto symGroupIt = symmetryGroupNodes.find(groupName);
        if (symGroupIt != symmetryGroupNodes.end()) {
            auto asfTree = symGroupIt->second->getASFTree();
            if (asfTree) {
                // Enforce symmetry constraints
                asfTree->enforceSymmetryConstraints();
            }
        }
    }
    
    // Clear affected cache
    clearAffectedCache();
    
    packed = true;
    
    return !hasOverlap();
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