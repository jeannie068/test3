/**
 * HBStarTree_Base.cpp
 * 
 * Implementation of the core structure and basic operations of the HBStarTree class
 * for analog placement with symmetry constraints.
 */

#include "HBStarTree.hpp"
#include <algorithm>
#include <random>
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
 * Clear the affected cache
 */
void HBStarTree::clearAffectedCache() {
    affectedModules.clear();
    affectedNodes.clear();
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
 * Gets all module nodes in the tree
 */
const map<string, shared_ptr<HBStarTreeNode>>& HBStarTree::getModuleNodes() const {
    return moduleNodes;
}

/**
 * Gets a random node from the tree
 */
shared_ptr<HBStarTreeNode> HBStarTree::randomNode(std::mt19937& rng) {
    auto& nodes = moduleNodes;
    if (nodes.empty()) return nullptr;
    
    std::uniform_int_distribution<int> dist(0, nodes.size() - 1);
    int index = dist(rng);
    auto it = nodes.begin();
    std::advance(it, index);
    return it->second;
}

/**
 * New function to implement global placement strategy with offset management
 */
void HBStarTree::improveGlobalPlacement() {
    // Phase 1: Place large modules first
    placeLargeModules();
    
    // Phase 2: Sort symmetry groups by area for better placement
    sortSymmetryGroupsByArea();
    
    // Now place symmetry groups in a more compact grid layout
    int maxRowWidth = MAX_ROW_WIDTH - 100; // Reduce max row width to ensure compactness
    int currentX = 0;
    int currentY = 0;
    int rowHeight = 0;
    int largeModulesMaxY = 0;
    
    // Find maximum Y for large modules
    for (const auto& pair : modules) {
        const auto& module = pair.second;
        
        if (!isInSymmetryGroup(pair.first) && 
            module->getWidth() * module->getHeight() > LARGE_MODULE_THRESHOLD) {
            largeModulesMaxY = std::max(largeModulesMaxY, module->getY() + module->getHeight());
        }
    }
    
    // Add a smaller buffer after large modules
    currentY = largeModulesMaxY + BUFFER;
    
    // Place symmetry groups in a grid to minimize overall width
    for (const auto& group : symmetryGroups) {
        // Get dimensions of this group
        auto dimensions = getSymmetryIslandDimensions(group);
        int width = dimensions.first;
        int height = dimensions.second;
        
        // Check if adding this group to current row would exceed width limit
        if (currentX + width > maxRowWidth && currentX > 0) {
            // Start a new row
            currentX = 0;
            currentY += rowHeight + BUFFER;
            rowHeight = 0;
        }
        
        // Get the hierarchy node for this group
        auto it = symmetryGroupNodes.find(group->getName());
        if (it == symmetryGroupNodes.end() || !it->second) {
            std::cerr << "Could not find hierarchy node for group: " << group->getName() << std::endl;
            continue;
        }
        
        // Position this symmetry group with exact coordinates
        positionSymmetryIsland(it->second, currentX, currentY);
        
        // Update for next placement
        currentX += width + BUFFER; // Reduced buffer
        rowHeight = std::max(rowHeight, height);
        
        std::cout << "Placed symmetry group " << group->getName() << " at (" 
                  << currentX - width - BUFFER << "," << currentY << ")" << std::endl;
    }
    
    // Phase 3: Place remaining regular modules in a compact grid
    currentY += rowHeight + BUFFER;
    currentX = 0;
    rowHeight = 0;
    
    // Collect regular modules (not in symmetry groups)
    std::vector<std::pair<std::string, std::shared_ptr<Module>>> regularModules;
    
    for (const auto& pair : modules) {
        const auto& moduleName = pair.first;
        const auto& module = pair.second;
        
        // Skip modules in symmetry groups
        if (isInSymmetryGroup(moduleName)) {
            continue;
        }
        
        // Skip large modules (already placed)
        if (module->getWidth() * module->getHeight() > LARGE_MODULE_THRESHOLD) {
            continue;
        }
        
        regularModules.push_back(pair);
    }
    
    // Sort regular modules by area (descending) for better packing
    std::sort(regularModules.begin(), regularModules.end(), 
        [](const auto& a, const auto& b) {
            return (a.second->getWidth() * a.second->getHeight()) > 
                   (b.second->getWidth() * b.second->getHeight());
        });
    
    // Place remaining modules in a grid with minimal gaps
    for (const auto& pair : regularModules) {
        const auto& module = pair.second;
        
        // Check if adding this module to current row would exceed width limit
        if (currentX + module->getWidth() > maxRowWidth && currentX > 0) {
            // Start a new row
            currentX = 0;
            currentY += rowHeight + BUFFER;
            rowHeight = 0;
        }
        
        // Position this module
        module->setPosition(currentX, currentY);
        
        // Update for next placement
        currentX += module->getWidth() + BUFFER; // Smaller buffer
        rowHeight = std::max(rowHeight, module->getHeight());
        
        std::cout << "Placed regular module " << pair.first << " at (" 
                  << module->getX() << "," << module->getY() << ")" << std::endl;
    }
    
    // Phase 4: Apply compaction to reduce any remaining space
    // compactPlacement();
    
    // Phase 5: Fix any overlaps that might have been introduced
    if (hasOverlap()) {
        std::cout << "Fixing any remaining overlaps after placement optimization..." << std::endl;
        fixGlobalOverlaps();
    }
    
    // Final enforcement of symmetry constraints
    for (const auto& symGroupPair : symmetryGroupNodes) {
        auto asfTree = symGroupPair.second->getASFTree();
        if (asfTree) {
            asfTree->enforceSymmetryConstraints();
        }
    }
    
    std::cout << "Global placement strategy completed" << std::endl;
}


/**
 * New function to sort symmetry groups by area
 */
void HBStarTree::sortSymmetryGroupsByArea() {
    // Sort symmetry groups by their area (largest first)
    std::sort(symmetryGroups.begin(), symmetryGroups.end(), 
        [this](const std::shared_ptr<SymmetryGroup>& a, const std::shared_ptr<SymmetryGroup>& b) {
            return getGroupArea(a->getName()) > getGroupArea(b->getName());
        });
}

/**
 * Position a symmetry island at the specified coordinates
 */
bool HBStarTree::positionSymmetryIsland(std::shared_ptr<HBStarTreeNode> hierarchyNode, int x, int y) {
    if (!hierarchyNode || hierarchyNode->getType() != HBNodeType::HIERARCHY) return false;
    
    auto asfTree = hierarchyNode->getASFTree();
    if (!asfTree) return false;
    
    // Find the bounding rectangle of the symmetry island
    int minX = std::numeric_limits<int>::max();
    int minY = std::numeric_limits<int>::max();
    int maxX = 0;
    int maxY = 0;
    
    for (const auto& pair : asfTree->getModules()) {
        const auto& module = pair.second;
        
        minX = std::min(minX, module->getX());
        minY = std::min(minY, module->getY());
        maxX = std::max(maxX, module->getX() + module->getWidth());
        maxY = std::max(maxY, module->getY() + module->getHeight());
    }
    
    // Calculate the offsets to shift the island
    int shiftX = x - minX;
    int shiftY = y - minY;
    
    // Update symmetry axis position
    if (asfTree->getSymmetryGroup()->getType() == SymmetryType::VERTICAL) {
        double oldAxis = asfTree->getSymmetryAxisPosition();
        asfTree->getSymmetryGroup()->setAxisPosition(oldAxis + shiftX);
    } else { // HORIZONTAL
        double oldAxis = asfTree->getSymmetryAxisPosition();
        asfTree->getSymmetryGroup()->setAxisPosition(oldAxis + shiftY);
    }
    
    // Shift all modules in the symmetry island
    for (const auto& pair : asfTree->getModules()) {
        auto module = pair.second;
        module->setPosition(module->getX() + shiftX, module->getY() + shiftY);
    }
    
    std::cout << "Positioned symmetry island " << hierarchyNode->getName() 
              << " at (" << x << "," << y << ") with shift (" 
              << shiftX << "," << shiftY << ")" << std::endl;
    
    // Update horizontal contour
    auto segments = asfTree->getContours().first->getSegments();
    for (const auto& segment : segments) {
        horizontalContour->addSegment(
            segment.start + shiftX,
            segment.end + shiftX,
            segment.height + shiftY
        );
    }
    
    // Update vertical contour
    auto vSegments = asfTree->getContours().second->getSegments();
    for (const auto& segment : vSegments) {
        verticalContour->addSegment(
            segment.start + shiftY,
            segment.end + shiftY,
            segment.height + shiftX
        );
    }
    
    return true;
}

/**
 * Get the height of a symmetry island
 */
int HBStarTree::getSymmetryIslandHeight(std::shared_ptr<SymmetryGroup> group) const {
    auto it = symmetryGroupNodes.find(group->getName());
    if (it == symmetryGroupNodes.end() || !it->second) return 0;
    
    auto asfTree = it->second->getASFTree();
    if (!asfTree) return 0;
    
    int minY = std::numeric_limits<int>::max();
    int maxY = 0;
    
    for (const auto& pair : asfTree->getModules()) {
        const auto& module = pair.second;
        minY = std::min(minY, module->getY());
        maxY = std::max(maxY, module->getY() + module->getHeight());
    }
    
    return (minY < maxY) ? (maxY - minY) : 0;
}

/**
 * Get the width of a symmetry island
 */
int HBStarTree::getSymmetryIslandWidth(std::shared_ptr<SymmetryGroup> group) const {
    auto it = symmetryGroupNodes.find(group->getName());
    if (it == symmetryGroupNodes.end() || !it->second) return 0;
    
    auto asfTree = it->second->getASFTree();
    if (!asfTree) return 0;
    
    int minX = std::numeric_limits<int>::max();
    int maxX = 0;
    
    for (const auto& pair : asfTree->getModules()) {
        const auto& module = pair.second;
        minX = std::min(minX, module->getX());
        maxX = std::max(maxX, module->getX() + module->getWidth());
    }
    
    return (minX < maxX) ? (maxX - minX) : 0;
}

/**
 * Get the total area of a symmetry group
 */
int HBStarTree::getGroupArea(const std::string& groupName) const {
    auto it = symmetryGroupNodes.find(groupName);
    if (it == symmetryGroupNodes.end() || !it->second) return 0;
    
    auto asfTree = it->second->getASFTree();
    if (!asfTree) return 0;
    
    int width = getSymmetryIslandWidth(asfTree->getSymmetryGroup());
    int height = getSymmetryIslandHeight(asfTree->getSymmetryGroup());
    
    return width * height;
}

/**
 * Shift an entire symmetry group by the specified amounts
 */
void HBStarTree::shiftSymmetryGroup(const std::string& groupName, int shiftX, int shiftY) {
    auto it = symmetryGroupNodes.find(groupName);
    if (it == symmetryGroupNodes.end() || !it->second) return;
    
    auto hierarchyNode = it->second;
    auto asfTree = hierarchyNode->getASFTree();
    
    if (!asfTree) return;
    
    // Update symmetry axis position
    if (asfTree->getSymmetryGroup()->getType() == SymmetryType::VERTICAL) {
        double oldAxis = asfTree->getSymmetryAxisPosition();
        asfTree->getSymmetryGroup()->setAxisPosition(oldAxis + shiftX);
    } else { // HORIZONTAL
        double oldAxis = asfTree->getSymmetryAxisPosition();
        asfTree->getSymmetryGroup()->setAxisPosition(oldAxis + shiftY);
    }
    
    // Shift all modules in the symmetry island
    for (const auto& pair : asfTree->getModules()) {
        auto module = pair.second;
        module->setPosition(module->getX() + shiftX, module->getY() + shiftY);
    }
    
    std::cout << "Shifted symmetry group " << groupName 
              << " by (" << shiftX << "," << shiftY << ")" << std::endl;
}

/**
 * Find which symmetry group a module belongs to
 */
std::string HBStarTree::findSymmetryGroupForModule(const std::string& moduleName) const {
    for (const auto& group : symmetryGroups) {
        // Check if module is in symmetry pairs
        for (const auto& pair : group->getSymmetryPairs()) {
            if (pair.first == moduleName || pair.second == moduleName) {
                return group->getName();
            }
        }
        
        // Check if module is self-symmetric
        for (const auto& name : group->getSelfSymmetric()) {
            if (name == moduleName) {
                return group->getName();
            }
        }
    }
    
    // Not found in any symmetry group
    return "";
}

/**
 * Find all overlapping module pairs in the design
 */
std::vector<std::pair<std::shared_ptr<Module>, std::shared_ptr<Module>>> 
HBStarTree::findAllOverlappingModulePairs() const {
    std::vector<std::pair<std::shared_ptr<Module>, std::shared_ptr<Module>>> overlappingPairs;
    
    // Collect all modules from all symmetry groups and regular modules
    std::vector<std::shared_ptr<Module>> allModules;
    
    // Add regular modules
    for (const auto& pair : modules) {
        allModules.push_back(pair.second);
    }
    
    // Add modules from symmetry groups
    for (const auto& pair : symmetryGroupNodes) {
        auto asfTree = pair.second->getASFTree();
        if (asfTree) {
            for (const auto& modPair : asfTree->getModules()) {
                allModules.push_back(modPair.second);
            }
        }
    }
    
    // Check all module pairs for overlaps
    for (size_t i = 0; i < allModules.size(); ++i) {
        auto& module1 = allModules[i];
        
        for (size_t j = i + 1; j < allModules.size(); ++j) {
            auto& module2 = allModules[j];
            
            // Skip if modules have the same name (duplicates)
            if (module1->getName() == module2->getName()) continue;
            
            if (module1->overlaps(*module2)) {
                overlappingPairs.push_back({module1, module2});
                
                std::cout << "Overlap detected between modules: " 
                          << module1->getName() << " (" 
                          << module1->getX() << "," << module1->getY() << "," 
                          << module1->getWidth() << "," << module1->getHeight() << ") and " 
                          << module2->getName() << " (" 
                          << module2->getX() << "," << module2->getY() << "," 
                          << module2->getWidth() << "," << module2->getHeight() << ")" << std::endl;
            }
        }
    }
    
    return overlappingPairs;
}

/**
 * Fix global overlaps between symmetry groups
 */
bool HBStarTree::fixGlobalOverlaps() {
    bool anyOverlapFixed = false;
    int fixIterations = 0;
    
    // Continue until no more overlaps or max iterations reached
    while (fixIterations < MAX_FIX_ITERATIONS) {
        auto overlappingPairs = findAllOverlaps();
        if (overlappingPairs.empty()) {
            std::cout << "No overlaps found, exiting fix iterations" << std::endl;
            break;
        }
        
        fixIterations++;
        bool fixedAny = false;
        
        std::cout << "Fixing overlaps: iteration " << fixIterations 
                  << " with " << overlappingPairs.size() << " overlapping pairs" << std::endl;
        
        // Sort overlapping pairs by overlap area (largest first)
        std::sort(overlappingPairs.begin(), overlappingPairs.end(),
            [this](const auto& a, const auto& b) {
                int areaA = getOverlapWidth(a.first, a.second) * getOverlapHeight(a.first, a.second);
                int areaB = getOverlapWidth(b.first, b.second) * getOverlapHeight(b.first, b.second);
                return areaA > areaB;
            });
        
        // Group overlaps by symmetry group
        std::map<std::string, std::vector<std::pair<std::shared_ptr<Module>, std::shared_ptr<Module>>>> groupedOverlaps;
        
        for (const auto& pair : overlappingPairs) {
            std::string group1 = findSymmetryGroupForModule(pair.first->getName());
            std::string group2 = findSymmetryGroupForModule(pair.second->getName());
            
            if (!group1.empty() && group1 == group2) {
                groupedOverlaps[group1].push_back(pair);
            }
        }
        
        // First fix overlaps within symmetry groups
        for (auto& entry : groupedOverlaps) {
            const std::string& groupName = entry.first;
            auto& pairs = entry.second;
            bool groupFixed = false;
            
            if (pairs.empty()) continue;
            
            // Attempt to resolve intra-group overlaps by adjusting the symmetry axis
            auto symGroupNode = getSymmetryGroupNode(groupName);
            if (symGroupNode) {
                auto asfTree = symGroupNode->getASFTree();
                if (asfTree) {
                    // Try adjusting symmetry axis position
                    double oldAxis = asfTree->getSymmetryAxisPosition();
                    // For vertical symmetry, try increasing the axis position
                    if (asfTree->getSymmetryGroup()->getType() == SymmetryType::VERTICAL) {
                        double newAxis = oldAxis + 5.0;  // Increase axis by 5 units
                        asfTree->getSymmetryGroup()->setAxisPosition(newAxis);
                        
                        // Recompute positions with new axis
                        if (asfTree->pack()) {
                            groupFixed = true;
                            fixedAny = true;
                        } else {
                            // Revert if unsuccessful
                            asfTree->getSymmetryGroup()->setAxisPosition(oldAxis);
                        }
                    } else {
                        // Similar logic for horizontal symmetry
                        double newAxis = oldAxis + 5.0;  // Increase axis by 5 units
                        asfTree->getSymmetryGroup()->setAxisPosition(newAxis);
                        
                        // Recompute positions with new axis
                        if (asfTree->pack()) {
                            groupFixed = true;
                            fixedAny = true;
                        } else {
                            // Revert if unsuccessful
                            asfTree->getSymmetryGroup()->setAxisPosition(oldAxis);
                        }
                    }
                    
                    // If still not fixed, try emergency recovery
                    if (!groupFixed) {
                        asfTree->emergencyRecovery();
                        if (!asfTree->checkForOverlaps()) {
                            groupFixed = true;
                            fixedAny = true;
                        }
                    }
                }
            }
            
            // If still not fixed, try shifting the entire group
            if (!groupFixed) {
                shiftSymmetryGroup(groupName, 0, EXTRA_LARGE_BUFFER);
                fixedAny = true;
            }
        }
        
        // Process remaining regular module overlaps
        for (const auto& pair : overlappingPairs) {
            std::string group1 = findSymmetryGroupForModule(pair.first->getName());
            std::string group2 = findSymmetryGroupForModule(pair.second->getName());
            
            // Skip if already handled by symmetry group fixes
            if (!group1.empty() && group1 == group2) {
                continue;
            }
            
            // CASE 1: Both are regular modules
            if (group1.empty() && group2.empty()) {
                // Determine which module to move (smaller one)
                auto moduleToMove = (pair.first->getArea() <= pair.second->getArea()) ? pair.first : pair.second;
                auto otherModule = (moduleToMove == pair.first) ? pair.second : pair.first;
                
                // Calculate overlap dimensions
                int overlapWidth = getOverlapWidth(pair.first, pair.second);
                int overlapHeight = getOverlapHeight(pair.first, pair.second);
                
                // Choose optimal direction to shift (minimum displacement)
                int shiftRight = 0, shiftDown = 0;
                
                // Calculate distances to each edge
                int distToLeft = moduleToMove->getX() - otherModule->getX();
                int distToRight = (otherModule->getX() + otherModule->getWidth()) - moduleToMove->getX();
                int distToTop = (otherModule->getY() + otherModule->getHeight()) - moduleToMove->getY();
                int distToBottom = moduleToMove->getY() - otherModule->getY();
                
                // Find minimum shift required
                int minDist = std::min({ 
                    abs(distToLeft - moduleToMove->getWidth()), 
                    abs(distToRight), 
                    abs(distToTop), 
                    abs(distToBottom - moduleToMove->getHeight()) 
                });
                
                // Determine optimal shift
                if (minDist == abs(distToLeft - moduleToMove->getWidth())) {
                    shiftRight = -(overlapWidth + 1);
                } else if (minDist == abs(distToRight)) {
                    shiftRight = overlapWidth + 1;
                } else if (minDist == abs(distToTop)) {
                    shiftDown = overlapHeight + 1;
                } else {
                    shiftDown = -(overlapHeight + 1);
                }
                
                // Apply the shift
                moduleToMove->setPosition(
                    moduleToMove->getX() + shiftRight,
                    moduleToMove->getY() + shiftDown
                );
                
                fixedAny = true;
            }
            // CASE 2: One module in symmetry group, one regular
            else if (!group1.empty() || !group2.empty()) {
                std::shared_ptr<Module> regularModule;
                std::string symGroupName;
                
                if (group1.empty()) {
                    regularModule = pair.first;
                    symGroupName = group2;
                } else {
                    regularModule = pair.second;
                    symGroupName = group1;
                }
                
                // Move the regular module away from the symmetry group
                int overlapHeight = getOverlapHeight(pair.first, pair.second);
                shiftModuleDown(regularModule, overlapHeight + BUFFER);
                
                fixedAny = true;
            }
            // CASE 3: Modules in different symmetry groups
            else if (!group1.empty() && !group2.empty() && group1 != group2) {
                // Determine which group to move (smaller area)
                std::string groupToMove = (getGroupArea(group1) <= getGroupArea(group2)) ? group1 : group2;
                
                // Shift the group vertically
                int overlapHeight = getOverlapHeight(pair.first, pair.second);
                shiftSymmetryGroup(groupToMove, 0, overlapHeight + LARGE_BUFFER);
                
                fixedAny = true;
            }
        }
        
        // If we couldn't fix any overlaps this iteration, try more drastic measures
        if (!fixedAny && !overlappingPairs.empty()) {
            // Spread out all symmetry groups
            int offsetY = EXTRA_LARGE_BUFFER;
            
            for (size_t i = 1; i < symmetryGroups.size(); ++i) {
                shiftSymmetryGroup(symmetryGroups[i]->getName(), 0, offsetY);
                offsetY += EXTRA_LARGE_BUFFER;
            }
            
            fixedAny = true;
        }
        
        // If we still couldn't fix any overlaps, we're done
        if (!fixedAny) {
            std::cerr << "Warning: Could not fix any more overlaps" << std::endl;
            break;
        }
        
        anyOverlapFixed = true;
    }
    
    if (fixIterations >= MAX_FIX_ITERATIONS) {
        std::cerr << "Warning: Reached maximum fix iterations (" << MAX_FIX_ITERATIONS << ")" << std::endl;
        
        // Last resort: try emergency recovery for all symmetry groups
        for (const auto& group : symmetryGroups) {
            auto symGroupNode = getSymmetryGroupNode(group->getName());
            if (symGroupNode) {
                auto asfTree = symGroupNode->getASFTree();
                if (asfTree) {
                    asfTree->emergencyRecovery();
                }
            }
        }
    }
    
    return anyOverlapFixed;
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

void HBStarTree::placeLargeModules() {
    std::vector<std::pair<std::string, std::shared_ptr<Module>>> largeModules;
    
    // Identify large modules from non-symmetry modules
    for (const auto& pair : modules) {
        const auto& moduleName = pair.first;
        const auto& module = pair.second;
        
        // Skip modules in symmetry groups
        if (isInSymmetryGroup(moduleName)) {
            continue;
        }
        
        // Check if this is a large module
        if (module->getWidth() * module->getHeight() > LARGE_MODULE_THRESHOLD) {
            largeModules.push_back(pair);
            std::cout << "Large module identified: " << moduleName << " with area " 
                      << (module->getWidth() * module->getHeight()) << std::endl;
        }
    }
    
    // Sort large modules by area (descending)
    std::sort(largeModules.begin(), largeModules.end(), 
        [](const auto& a, const auto& b) {
            return (a.second->getWidth() * a.second->getHeight()) > 
                   (b.second->getWidth() * b.second->getHeight());
        });
    
    // Position large modules at the origin with extra spacing between them
    int xPosition = 0;
    int yPosition = 0;
    int rowHeight = 0;
    
    for (const auto& pair : largeModules) {
        const auto& module = pair.second;
        
        // Check if adding this module would exceed the row width
        if (xPosition + module->getWidth() > MAX_ROW_WIDTH && xPosition > 0) {
            // Start a new row
            xPosition = 0;
            yPosition += rowHeight + EXTRA_LARGE_BUFFER;
            rowHeight = 0;
        }
        
        // Position this module
        module->setPosition(xPosition, yPosition);
        
        // Update for next placement
        xPosition += module->getWidth() + EXTRA_LARGE_BUFFER;
        rowHeight = std::max(rowHeight, module->getHeight());
        
        std::cout << "Placed large module " << pair.first << " at (" 
                  << module->getX() << "," << module->getY() << ")" << std::endl;
    }
}

/**
 * Find all overlapping module pairs in the design with enhanced debugging
 */
std::vector<std::pair<std::shared_ptr<Module>, std::shared_ptr<Module>>> 
HBStarTree::findAllOverlaps() const {
    std::vector<std::pair<std::shared_ptr<Module>, std::shared_ptr<Module>>> overlappingPairs;
    
    // Collect all modules from all symmetry groups and regular modules
    std::vector<std::shared_ptr<Module>> allModules;
    
    // Add regular modules
    for (const auto& pair : modules) {
        allModules.push_back(pair.second);
    }
    
    // Add modules from symmetry groups
    for (const auto& pair : symmetryGroupNodes) {
        auto asfTree = pair.second->getASFTree();
        if (asfTree) {
            for (const auto& modPair : asfTree->getModules()) {
                allModules.push_back(modPair.second);
            }
        }
    }
    
    // Check all module pairs for overlaps
    for (size_t i = 0; i < allModules.size(); ++i) {
        auto& module1 = allModules[i];
        
        for (size_t j = i + 1; j < allModules.size(); ++j) {
            auto& module2 = allModules[j];
            
            // Skip if modules have the same name (duplicates)
            if (module1->getName() == module2->getName()) continue;
            
            // Check for overlap
            if (module1->overlaps(*module2)) {
                overlappingPairs.push_back({module1, module2});
                
                // Calculate overlap areas for detailed debugging
                int overlapWidth = getOverlapWidth(module1, module2);
                int overlapHeight = getOverlapHeight(module1, module2);
                int overlapArea = overlapWidth * overlapHeight;
                
                std::cout << "CRITICAL OVERLAP: " << module1->getName() << " and " 
                          << module2->getName() << " overlap by " << overlapArea 
                          << " square units" << std::endl;
                std::cout << "  " << module1->getName() << ": (" << module1->getX() 
                          << "," << module1->getY() << ") size=" << module1->getWidth() 
                          << "x" << module1->getHeight() << std::endl;
                std::cout << "  " << module2->getName() << ": (" << module2->getX() 
                          << "," << module2->getY() << ") size=" << module2->getWidth() 
                          << "x" << module2->getHeight() << std::endl;
            }
        }
    }
    
    std::cout << "Found " << overlappingPairs.size() << " overlapping module pairs" << std::endl;
    return overlappingPairs;
}

/**
 * Check if a module is a regular module (not in a symmetry group)
 */
bool HBStarTree::isRegularModule(const std::string& moduleName) const {
    return !isInSymmetryGroup(moduleName);
}

/**
 * Check if a module is in any symmetry group
 */
bool HBStarTree::isInSymmetryGroup(const std::string& moduleName) const {
    for (const auto& group : symmetryGroups) {
        // Check symmetry pairs
        for (const auto& pair : group->getSymmetryPairs()) {
            if (pair.first == moduleName || pair.second == moduleName) {
                return true;
            }
        }
        
        // Check self-symmetric modules
        for (const auto& name : group->getSelfSymmetric()) {
            if (name == moduleName) {
                return true;
            }
        }
    }
    
    return false;
}

/**
 * Get the area of a module
 */
int HBStarTree::getModuleArea(const std::string& moduleName) const {
    auto it = modules.find(moduleName);
    if (it == modules.end()) return 0;
    
    auto module = it->second;
    return module->getWidth() * module->getHeight();
}

/**
 * Get the width of the overlap between two modules
 */
int HBStarTree::getOverlapWidth(const std::shared_ptr<Module>& mod1, 
                              const std::shared_ptr<Module>& mod2) const {
    int left = std::max(mod1->getX(), mod2->getX());
    int right = std::min(mod1->getX() + mod1->getWidth(), mod2->getX() + mod2->getWidth());
    
    return (right > left) ? (right - left) : 0;
}

/**
 * Get the height of the overlap between two modules
 */
int HBStarTree::getOverlapHeight(const std::shared_ptr<Module>& mod1, 
                               const std::shared_ptr<Module>& mod2) const {
    int bottom = std::max(mod1->getY(), mod2->getY());
    int top = std::min(mod1->getY() + mod1->getHeight(), mod2->getY() + mod2->getHeight());
    
    return (top > bottom) ? (top - bottom) : 0;
}

/**
 * Shift a module horizontally
 */
void HBStarTree::shiftModuleRight(const std::shared_ptr<Module>& module, int shift) {
    if (!module) return;
    
    int newX = module->getX() + shift;
    // Ensure we don't go below 0
    newX = std::max(0, newX);
    
    std::cout << "Shifting module " << module->getName() << " horizontally from " 
              << module->getX() << " to " << newX << std::endl;
    
    module->setPosition(newX, module->getY());
}

/**
 * Shift a module vertically
 */
void HBStarTree::shiftModuleDown(const std::shared_ptr<Module>& module, int shift) {
    if (!module) return;
    
    int newY = module->getY() + shift;
    // Ensure we don't go below 0
    newY = std::max(0, newY);
    
    std::cout << "Shifting module " << module->getName() << " vertically from " 
              << module->getY() << " to " << newY << std::endl;
    
    module->setPosition(module->getX(), newY);
}

/**
 * Get the dimensions of a symmetry island
 */
std::pair<int, int> HBStarTree::getSymmetryIslandDimensions(
    const std::shared_ptr<SymmetryGroup>& group) const {
    auto it = symmetryGroupNodes.find(group->getName());
    if (it == symmetryGroupNodes.end() || !it->second) {
        return {0, 0};
    }
    
    auto asfTree = it->second->getASFTree();
    if (!asfTree) {
        return {0, 0};
    }
    
    int minX = std::numeric_limits<int>::max();
    int minY = std::numeric_limits<int>::max();
    int maxX = 0;
    int maxY = 0;
    
    for (const auto& pair : asfTree->getModules()) {
        const auto& module = pair.second;
        
        minX = std::min(minX, module->getX());
        minY = std::min(minY, module->getY());
        maxX = std::max(maxX, module->getX() + module->getWidth());
        maxY = std::max(maxY, module->getY() + module->getHeight());
    }
    
    // Return width and height
    return {maxX - minX, maxY - minY};
}

/**
 * Verify that no overlaps exist in the final placement
 */
bool HBStarTree::verifyNoOverlaps() const {
    auto overlaps = findAllOverlaps();
    return overlaps.empty();
}