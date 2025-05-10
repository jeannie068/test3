/**
 * HBStarTree.hpp
 * 
 * This file defines the HBStarTree class, which is a hierarchical framework
 * that can simultaneously optimize the placement with both symmetry islands
 * and non-symmetric modules.
 * 
 * Based on the paper "Analog Placement Based on Symmetry-Island Formulation".
 */

#pragma once

#include <memory>
#include <vector>
#include <map>
#include <string>
#include <utility>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <random>

#include "HBStarTreeNode.hpp"
#include "ASFBStarTree.hpp"
#include "Module.hpp"
#include "SymmetryConstraint.hpp"
#include "../utils/Contour.hpp"
using namespace std;

// Forward declaration to avoid circular dependency
class ASFBStarTree;

class HBStarTree {
private:
    shared_ptr<HBStarTreeNode> root;
    
    map<string, shared_ptr<Module>> modules; // All modules
    vector<shared_ptr<SymmetryGroup>> symmetryGroups; // All symmetry groups 
    
    // Map from symmetry group name to its hierarchy node
    map<string, shared_ptr<HBStarTreeNode>> symmetryGroupNodes;
    
    // Map from module name to its node in the tree
    map<string, shared_ptr<HBStarTreeNode>> moduleNodes;
    
    // Contour for efficient packing
    shared_ptr<Contour> horizontalContour;
    shared_ptr<Contour> verticalContour;
    
    // Total area of the placement
    int totalArea;
    
    // Flag indicating whether the tree has been packed
    bool packed;
    
    // Cache of affected modules and nodes for incremental packing
    unordered_set<string> affectedModules;
    unordered_set<shared_ptr<HBStarTreeNode>> affectedNodes;

    // New functions for global placement management
    void improveGlobalPlacement();
    bool positionSymmetryIsland(std::shared_ptr<HBStarTreeNode> hierarchyNode, int x, int y);
    int getSymmetryIslandHeight(std::shared_ptr<SymmetryGroup> group) const;
    int getSymmetryIslandWidth(std::shared_ptr<SymmetryGroup> group) const;
    void shiftSymmetryGroup(const std::string& groupName, int shiftX, int shiftY);
    bool fixGlobalOverlaps();
    void sortSymmetryGroupsByArea();
    std::vector<std::pair<std::shared_ptr<Module>, std::shared_ptr<Module>>> findAllOverlappingModulePairs() const;
    std::string findSymmetryGroupForModule(const std::string& moduleName) const;
    int getGroupArea(const std::string& groupName) const;
    
    // Internal helper methods
    void updateContourNodes();
    shared_ptr<HBStarTreeNode> findNearestContourNode(shared_ptr<HBStarTreeNode> node) const;
    shared_ptr<HBStarTreeNode> findLeftmostSkewedChild(shared_ptr<HBStarTreeNode> node) const;
    void constructSymmetryIslands();
    void constructInitialTreeStructure();
    void clearTree();
    
    // Helper methods for incremental packing
    bool packSubtree(shared_ptr<HBStarTreeNode> node);
    bool packNode(shared_ptr<HBStarTreeNode> node);
    bool packSymmetryIsland(shared_ptr<HBStarTreeNode> hierarchyNode);
    void markAffectedSubtree(shared_ptr<HBStarTreeNode> node);
    void clearAffectedCache();

    // Constants for placement parameters
    static constexpr int BUFFER = 1;
    static constexpr int LARGE_BUFFER = 2;
    static constexpr int EXTRA_LARGE_BUFFER = 5;
    static constexpr int MAX_ROW_WIDTH = 2000;
    static constexpr int MAX_FIX_ITERATIONS = 15;
    static constexpr int LARGE_MODULE_THRESHOLD = 10000; // Area threshold for large modules
    
    // New helper functions
    bool isRegularModule(const std::string& moduleName) const;
    bool isInSymmetryGroup(const std::string& moduleName) const;
    int getModuleArea(const std::string& moduleName) const;
    int getOverlapWidth(const std::shared_ptr<Module>& mod1, const std::shared_ptr<Module>& mod2) const;
    int getOverlapHeight(const std::shared_ptr<Module>& mod1, const std::shared_ptr<Module>& mod2) const;
    void shiftModuleRight(const std::shared_ptr<Module>& module, int shift);
    void shiftModuleDown(const std::shared_ptr<Module>& module, int shift);
    std::pair<int, int> getSymmetryIslandDimensions(const std::shared_ptr<SymmetryGroup>& group) const;
    void placeLargeModules();
    bool verifyNoOverlaps() const;
    std::vector<std::pair<std::shared_ptr<Module>, std::shared_ptr<Module>>> findAllOverlaps() const;
    
public:
    /**
     * Constructor
     */
    HBStarTree();
    ~HBStarTree();
    
    /* Adds a module to the tree */
    void addModule(shared_ptr<Module> module);
    
    /* Adds a symmetry group to the tree */
    void addSymmetryGroup(shared_ptr<SymmetryGroup> group);
    
    /* Constructs an initial HB*-tree */
    void constructInitialTree();
    
    /**
     * Calculates the coordinates of all modules by packing the HB*-tree
     * 
     * @return True if packing was successful, false otherwise
     */
    bool pack();
    
    /**
     * Incrementally updates the packing after a perturbation
     * 
     * @return True if incremental packing was successful, false otherwise
     */
    bool incrementalPack();
    
    /**
     * Returns whether the tree has been packed
     */
    bool isPacked() const;
    
    /**
     * Performs a rotation operation on a module
     * 
     * @param moduleName Name of the module to rotate
     * @return True if the rotation was successful, false otherwise
     */
    bool rotateModule(const string& moduleName);
    
    /**
     * Moves a node to a new position in the tree
     * 
     * @param nodeName Name of the node to move
     * @param newParentName Name of the new parent node
     * @param asLeftChild True if the node should be the left child, false for right child
     * @return True if the move was successful, false otherwise
     */
    bool moveNode(const string& nodeName, 
                  const string& newParentName, 
                  bool asLeftChild);
    
    /**
     * Swaps two nodes in the tree
     * 
     * @param nodeName1 Name of the first node
     * @param nodeName2 Name of the second node
     * @return True if the swap was successful, false otherwise
     */
    bool swapNodes(const string& nodeName1, const string& nodeName2);
    
    /**
     * Changes the representative of a symmetry pair in a symmetry group
     * 
     * @param symmetryGroupName Name of the symmetry group
     * @param moduleName Name of the module in the symmetry pair
     * @return True if the change was successful, false otherwise
     */
    bool changeRepresentative(const string& symmetryGroupName, const string& moduleName);
    
    /**
     * Converts the symmetry type of a symmetry group
     * 
     * @param symmetryGroupName Name of the symmetry group
     * @return True if the conversion was successful, false otherwise
     */
    bool convertSymmetryType(const string& symmetryGroupName);

    /**
     * Checks if the placement has any module overlaps
     * 
     * @return True if there are overlaps, false otherwise
     */
    bool hasOverlap() const;
    
    /* Returns the total area of the placement */
    int getArea() const;
    
    /**
     * Returns the total wire length of the placement
     * 
     * @return Total wire length of the placement
     */
    int getWireLength() const;
    
    /**
     * Gets the root node of the HB*-tree
     * 
     * @return Root node of the HB*-tree
     */
    shared_ptr<HBStarTreeNode> getRoot() const;
    
    /**
     * Gets all modules in the design
     * 
     * @return Map of module names to modules
     */
    const map<string, shared_ptr<Module>>& getModules() const;
    
    /**
     * Gets all symmetry groups in the design
     * 
     * @return Vector of symmetry groups
     */
    const vector<shared_ptr<SymmetryGroup>>& getSymmetryGroups() const;
    
    /**
     * Gets the module node with the given name
     * 
     * @param moduleName Name of the module
     * @return Module node with the given name, or nullptr if not found
     */
    shared_ptr<HBStarTreeNode> getModuleNode(const string& moduleName) const;

    /**
     * Gets all module nodes in the tree
     * 
     * @return Map of module names to nodes
     */
    const map<string, shared_ptr<HBStarTreeNode>>& getModuleNodes() const;
    
    /**
     * Gets the symmetry group node with the given name
     * 
     * @param symmetryGroupName Name of the symmetry group
     * @return Symmetry group node with the given name, or nullptr if not found
     */
    shared_ptr<HBStarTreeNode> getSymmetryGroupNode(const string& symmetryGroupName) const;
    
    /**
     * Marks a module as affected, needing repacking
     * 
     * @param moduleName Name of the affected module
     */
    void markModuleAffected(const string& moduleName);
    
    /**
     * Marks a symmetry group as affected, needing repacking
     * 
     * @param symmetryGroupName Name of the affected symmetry group
     */
    void markSymmetryGroupAffected(const string& symmetryGroupName);
    
    /**
     * Marks a node as affected, needing repacking
     * 
     * @param node The affected node
     */
    void markNodeAffected(shared_ptr<HBStarTreeNode> node);

    /**
     * Gets a random node from the tree
     * 
     * @param rng Random number generator
     * @return Random node from the tree, or nullptr if empty
     */
    shared_ptr<HBStarTreeNode> randomNode(std::mt19937& rng);

    /**
     * NEW: Method to rebuild the global contours and update area calculation
     */
    void rebuildGlobalContours();
    
    /**
     * Creates a deep copy of this HB*-tree
     * 
     * @return A new HB*-tree that is a deep copy of this one
     */
    shared_ptr<HBStarTree> clone() const;
};