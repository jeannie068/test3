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
    // Phase 1: Place large modules with minimum spacing
    placeLargeModules();
    
    // Organize symmetry groups in decreasing order of area
    sortSymmetryGroupsByArea();
    
    // Create contour-based placement for optimal compaction
    std::shared_ptr<Contour> globalHorizontalContour = std::make_shared<Contour>();
    std::shared_ptr<Contour> globalVerticalContour = std::make_shared<Contour>();
    
    // Initialize global contours with reasonable width
    globalHorizontalContour->addSegment(0, 10000, 0);
    globalVerticalContour->addSegment(0, 10000, 0);
    
    // Find maximum Y of large modules to start placement
    int largeModulesMaxY = 0;
    for (const auto& pair : modules) {
        const auto& module = pair.second;
        
        if (!isInSymmetryGroup(pair.first) && 
            module->getWidth() * module->getHeight() > LARGE_MODULE_THRESHOLD) {
            largeModulesMaxY = std::max(largeModulesMaxY, module->getY() + module->getHeight());
        }
    }
    
    // Add minimal buffer after large modules
    int startY = largeModulesMaxY + 1;
    
    // Phase 2: Place symmetry groups using contour-based approach
    for (const auto& group : symmetryGroups) {
        // Get dimensions of this group
        auto dimensions = getSymmetryIslandDimensions(group);
        int width = dimensions.first;
        int height = dimensions.second;
        
        // Find best position along horizontal contour for minimal area
        int bestX = 0;
        int bestY = std::numeric_limits<int>::max();
        int searchStep = 5; // Step size for efficiency
        
        for (int testX = 0; testX < 5000; testX += searchStep) {
            int yHeight = globalHorizontalContour->getHeight(testX, testX + width);
            if (yHeight < bestY) {
                bestX = testX;
                bestY = yHeight;
            }
        }
        
        // Ensure we don't place below start Y
        bestY = std::max(bestY, startY);
        
        // Get the hierarchy node for this group
        auto it = symmetryGroupNodes.find(group->getName());
        if (it == symmetryGroupNodes.end() || !it->second) {
            std::cerr << "Could not find hierarchy node for group: " << group->getName() << std::endl;
            continue;
        }
        
        // Position this symmetry group with exact coordinates
        positionSymmetryIsland(it->second, bestX, bestY);
        
        // Update contours with this symmetry island
        globalHorizontalContour->addSegment(bestX, bestX + width, bestY + height);
        globalVerticalContour->addSegment(bestY, bestY + height, bestX + width);
        
        std::cout << "Placed symmetry group " << group->getName() << " at (" 
                  << bestX << "," << bestY << ")" << std::endl;
    }
    
    // Phase 3: Place remaining regular modules using contour-based approach
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
    
    // Place regular modules using contour-based approach
    for (const auto& pair : regularModules) {
        const auto& module = pair.second;
        int width = module->getWidth();
        int height = module->getHeight();
        
        // Find best position along horizontal contour
        int bestX = 0;
        int bestY = std::numeric_limits<int>::max();
        
        for (int testX = 0; testX < 5000; testX += 5) {
            int yHeight = globalHorizontalContour->getHeight(testX, testX + width);
            if (yHeight < bestY) {
                bestX = testX;
                bestY = yHeight;
            }
        }
        
        // Position module
        module->setPosition(bestX, bestY);
        
        // Update contours
        globalHorizontalContour->addSegment(bestX, bestX + width, bestY + height);
        globalVerticalContour->addSegment(bestY, bestY + height, bestX + width);
        
        std::cout << "Placed regular module " << pair.first << " at (" 
                  << module->getX() << "," << module->getY() << ")" << std::endl;
    }
    
    // Phase 4: Fix any overlaps that might have been introduced
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
    
    // Rebuild global contours for accurate area calculation
    rebuildGlobalContours();
    
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
 * FIX: Remove old contours before adding new ones
 */
// In HBStarTree::positionSymmetryIsland - improved removal of old contours
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
    
    // Get the island's original contours before shifting
    auto oldContours = asfTree->getContours();
    auto oldHorizontalSegments = oldContours.first->getSegments();
    auto oldVerticalSegments = oldContours.second->getSegments();
    
    // Properly remove old contour segments from global contours
    for (const auto& segment : oldHorizontalSegments) {
        // Create a temporary segment with the same parameters to avoid modifying the original
        ContourSegment tempSegment(segment.start, segment.end, segment.height);
        // Use removeSegment method to remove this segment
        horizontalContour->removeSegment(tempSegment);
    }
    
    for (const auto& segment : oldVerticalSegments) {
        ContourSegment tempSegment(segment.start, segment.end, segment.height);
        verticalContour->removeSegment(tempSegment);
    }
    
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
    
    // Rebuild the internal contours of the ASF-B*-tree
    asfTree->rebuildLocalContours();
    
    // Get the updated contours after shifting and rebuilding
    auto newContours = asfTree->getContours();
    auto segments = newContours.first->getSegments();
    auto vSegments = newContours.second->getSegments();
    
    // Update horizontal contour with the new segments
    for (const auto& segment : segments) {
        horizontalContour->addSegment(
            segment.start,
            segment.end,
            segment.height
        );
    }
    
    // Update vertical contour with the new segments
    for (const auto& segment : vSegments) {
        verticalContour->addSegment(
            segment.start,
            segment.end,
            segment.height
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
    const int MAX_FIX_ITERATIONS = 50; // Increased from current value
    
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
        
        // First fix overlaps within symmetry groups by enforcing constraints
        for (auto& entry : groupedOverlaps) {
            const std::string& groupName = entry.first;
            auto& pairs = entry.second;
            bool groupFixed = false;
            
            if (pairs.empty()) continue;
            
            // Fix intra-group overlaps by enforcing symmetry constraints
            auto symGroupNode = getSymmetryGroupNode(groupName);
            if (symGroupNode) {
                auto asfTree = symGroupNode->getASFTree();
                if (asfTree) {
                    // Apply symmetry constraint enforcement to fix overlaps
                    asfTree->enforceSymmetryConstraints();
                    
                    // Check if it resolved overlaps
                    if (!asfTree->checkForOverlaps()) {
                        groupFixed = true;
                        fixedAny = true;
                        std::cout << "Fixed overlaps in symmetry group " << groupName 
                                  << " by enforcing constraints" << std::endl;
                    } else {
                        // If not fixed, try vertical separation
                        for (size_t i = 1; i < pairs.size(); i++) {
                            auto prevModule = pairs[i-1].first;
                            auto currModule = pairs[i].first;
                            
                            // Separate vertically with minimal gap
                            int newY = prevModule->getY() + prevModule->getHeight() + 1;
                            currModule->setPosition(currModule->getX(), newY);
                            
                            // Also move the symmetric pair
                            auto symModule = pairs[i].second;
                            symModule->setPosition(symModule->getX(), newY);
                            
                            groupFixed = true;
                            fixedAny = true;
                        }
                    }
                    
                    // Final enforcement to ensure symmetry is maintained
                    asfTree->enforceSymmetryConstraints();
                }
            }
            
            // If that didn't work, shift the entire group
            if (!groupFixed) {
                shiftSymmetryGroup(groupName, 0, 10); // Minimal shift of 10 units vertically
                fixedAny = true;
            }
        }
        
        // Process remaining module overlaps with minimal displacements
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
                
                // Choose minimal displacement direction
                if (overlapWidth <= overlapHeight) {
                    // Horizontal shift with minimal displacement
                    int shiftX = 0;
                    
                    if (moduleToMove->getX() < otherModule->getX()) {
                        // Move left
                        shiftX = -(overlapWidth + 1);
                    } else {
                        // Move right
                        shiftX = overlapWidth + 1;
                    }
                    
                    moduleToMove->setPosition(
                        moduleToMove->getX() + shiftX,
                        moduleToMove->getY()
                    );
                } else {
                    // Vertical shift with minimal displacement
                    int shiftY = 0;
                    
                    if (moduleToMove->getY() < otherModule->getY()) {
                        // Move up
                        shiftY = -(overlapHeight + 1);
                    } else {
                        // Move down
                        shiftY = overlapHeight + 1;
                    }
                    
                    moduleToMove->setPosition(
                        moduleToMove->getX(),
                        moduleToMove->getY() + shiftY
                    );
                }
                
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
                
                // Move the regular module away with minimal displacement
                int overlapHeight = getOverlapHeight(pair.first, pair.second);
                // Add just 1 unit more than needed
                shiftModuleDown(regularModule, overlapHeight + 1);
                
                fixedAny = true;
            }
            // CASE 3: Modules in different symmetry groups
            else if (!group1.empty() && !group2.empty() && group1 != group2) {
                // Determine which group to move (smaller area)
                std::string groupToMove = (getGroupArea(group1) <= getGroupArea(group2)) ? group1 : group2;
                
                // Shift the group vertically with minimal displacement
                int overlapHeight = getOverlapHeight(pair.first, pair.second);
                shiftSymmetryGroup(groupToMove, 0, overlapHeight + 1);
                
                fixedAny = true;
            }
        }
        
        // If we couldn't fix any overlaps this iteration, try more drastic measures
        if (!fixedAny && !overlappingPairs.empty()) {
            // Spread out all symmetry groups with minimal spacing
            int offsetY = 10;
            
            for (size_t i = 1; i < symmetryGroups.size(); ++i) {
                shiftSymmetryGroup(symmetryGroups[i]->getName(), 0, offsetY);
                offsetY += 10;
            }
            
            fixedAny = true;
        }
        
        // If we still couldn't fix any overlaps, we're done
        if (!fixedAny) {
            std::cerr << "Warning: Could not fix any more overlaps" << std::endl;
            break;
        }
        
        anyOverlapFixed = true;
        
        // Check if all overlaps are fixed
        if (!hasOverlap()) {
            std::cout << "All overlaps successfully fixed in iteration " << fixIterations << std::endl;
            break;
        }
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
                    // Enforce constraints after recovery
                    asfTree->enforceSymmetryConstraints();
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
        
        // std::cout << "Placed large module " << pair.first << " at (" 
        //           << module->getX() << "," << module->getY() << ")" << std::endl;
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
                
                /*
                std::cout << "CRITICAL OVERLAP: " << module1->getName() << " and " 
                          << module2->getName() << " overlap by " << overlapArea 
                          << " square units" << std::endl;
                std::cout << "  " << module1->getName() << ": (" << module1->getX() 
                          << "," << module1->getY() << ") size=" << module1->getWidth() 
                          << "x" << module1->getHeight() << std::endl;
                std::cout << "  " << module2->getName() << ": (" << module2->getX() 
                          << "," << module2->getY() << ") size=" << module2->getWidth() 
                          << "x" << module2->getHeight() << std::endl;
                */
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
    
    // std::cout << "Shifting module " << module->getName() << " horizontally from " 
    //           << module->getX() << " to " << newX << std::endl;
    
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
    
    // std::cout << "Shifting module " << module->getName() << " vertically from " 
    //           << module->getY() << " to " << newY << std::endl;
    
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

/**
 * NEW: Method to rebuild the global contours and update area calculation
 */
void HBStarTree::rebuildGlobalContours() {
    // Clear global contours
    horizontalContour->clear();
    verticalContour->clear();
    
    // Find maximum dimensions for reasonable width/height limits
    int maxWidth = 0;
    int maxHeight = 0;
    
    // Examine all modules to find the layout dimensions
    for (const auto& pair : modules) {
        const auto& module = pair.second;
        maxWidth = std::max(maxWidth, module->getX() + module->getWidth());
        maxHeight = std::max(maxHeight, module->getY() + module->getHeight());
    }
    
    // Also check modules in symmetry groups
    for (const auto& pair : symmetryGroupNodes) {
        auto asfTree = pair.second->getASFTree();
        if (asfTree) {
            for (const auto& modPair : asfTree->getModules()) {
                const auto& module = modPair.second;
                maxWidth = std::max(maxWidth, module->getX() + module->getWidth());
                maxHeight = std::max(maxHeight, module->getY() + module->getHeight());
            }
        }
    }
    
    // Add margin for safety
    maxWidth = maxWidth * 2;
    maxHeight = maxHeight * 2;
    
    // Initialize horizontal contour with reasonable width
    horizontalContour->initializeSentinel(maxWidth + 100);
    
    // Initialize vertical contour
    verticalContour->initializeSentinel(maxHeight + 100);
    
    // Add all modules to the contours
    for (const auto& pair : modules) {
        const auto& module = pair.second;
        horizontalContour->addSegment(
            module->getX(), 
            module->getX() + module->getWidth(), 
            module->getY() + module->getHeight()
        );
        
        verticalContour->addSegment(
            module->getY(), 
            module->getY() + module->getHeight(), 
            module->getX() + module->getWidth()
        );
    }
    
    // Add all symmetry island contours
    for (const auto& pair : symmetryGroupNodes) {
        auto asfTree = pair.second->getASFTree();
        if (asfTree) {
            auto contours = asfTree->getContours();
            auto hSegs = contours.first->getSegments();
            auto vSegs = contours.second->getSegments();
            
            for (const auto& seg : hSegs) {
                horizontalContour->addSegment(seg.start, seg.end, seg.height);
            }
            
            for (const auto& seg : vSegs) {
                verticalContour->addSegment(seg.start, seg.end, seg.height);
            }
        }
    }
    
    // Recalculate total area
    int maxX = 0, maxY = 0;
    
    for (const auto& pair : modules) {
        const auto& module = pair.second;
        maxX = std::max(maxX, module->getX() + module->getWidth());
        maxY = std::max(maxY, module->getY() + module->getHeight());
    }
    
    for (const auto& pair : symmetryGroupNodes) {
        auto asfTree = pair.second->getASFTree();
        if (asfTree) {
            for (const auto& modPair : asfTree->getModules()) {
                const auto& module = modPair.second;
                maxX = std::max(maxX, module->getX() + module->getWidth());
                maxY = std::max(maxY, module->getY() + module->getHeight());
            }
        }
    }
    
    totalArea = maxX * maxY;
    std::cout << "Rebuilt global contours and updated area to: " << totalArea << std::endl;
}