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