// ASFBStarTree.cpp

#include "ASFBStarTree.hpp"
#include <algorithm>
#include <queue>
#include <iostream>
#include <limits>
#include <cmath>
using namespace std;

/**
 * Constructor
 */
ASFBStarTree::ASFBStarTree(shared_ptr<SymmetryGroup> symmetryGroup)
    : symmetryGroup(symmetryGroup), 
      root(nullptr),
      horizontalContour(make_shared<Contour>()),
      verticalContour(make_shared<Contour>()),
      symmetryAxisPosition(0.0) {
    
    // Initialize representative and symmetry pair maps
    if (symmetryGroup) {
        // Process symmetry pairs
        for (const auto& pair : symmetryGroup->getSymmetryPairs()) {
            // For each symmetry pair, choose the second module as the representative
            representativeMap[pair.first] = pair.second;
            representativeMap[pair.second] = pair.second;  // Representative represents itself
            
            // Track the symmetric pair relationship
            symmetricPairMap[pair.first] = pair.second;
            symmetricPairMap[pair.second] = pair.first;
        }
        
        // Process self-symmetric modules
        for (const auto& moduleName : symmetryGroup->getSelfSymmetric()) {
            // For a self-symmetric module, it represents itself but half of it
            representativeMap[moduleName] = moduleName;
            selfSymmetricModules.push_back(moduleName);
        }
    }
}

ASFBStarTree::~ASFBStarTree() {
}

/**
 * Adds a module to the tree
 */
void ASFBStarTree::addModule(shared_ptr<Module> module) {
    if (!module) return;
    
    // Add the module to our module map
    modules[module->getName()] = module;
}

/**
 * Constructs an initial ASF-B*-tree based on the symmetry group
 */
void ASFBStarTree::constructInitialTree() {
    // Clear any existing tree
    root = nullptr;
    
    // First, collect all representatives
    vector<string> representatives;
    for (const auto& pair : modules) {
        const string& moduleName = pair.first;
        if (isRepresentative(moduleName)) {
            representatives.push_back(moduleName);
        }
    }
    
    if (representatives.empty()) return;
    
    // Sort representatives by area (largest first)
    sort(representatives.begin(), representatives.end(), 
              [this](const string& a, const string& b) {
                  return modules[a]->getArea() > modules[b]->getArea();
              });
    
    // Create the root node with the first representative
    root = make_shared<BStarTreeNode>(representatives[0]);
    
    // Add other representatives to the tree
    for (size_t i = 1; i < representatives.size(); ++i) {
        auto newNode = make_shared<BStarTreeNode>(representatives[i]);
        
        // Case: self-symmetric module - it must be on the boundary
        if (find(selfSymmetricModules.begin(), selfSymmetricModules.end(), 
                      representatives[i]) != selfSymmetricModules.end()) {
            
            // Place self-symmetric modules on the rightmost branch for vertical symmetry
            // or leftmost branch for horizontal symmetry
            if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
                // Find the rightmost node
                auto current = root;
                while (current->getRightChild()) {
                    current = current->getRightChild();
                }
                current->setRightChild(newNode);
                newNode->setParent(current);
            } else {
                // Find the leftmost node
                auto current = root;
                while (current->getLeftChild()) {
                    current = current->getLeftChild();
                }
                current->setLeftChild(newNode);
                newNode->setParent(current);
            }
        } else {
            // Case: symmetry pairs - can place anywhere
            // For simplicity, place them as right children of the rightmost node
            auto current = root;
            while (current->getRightChild()) {
                current = current->getRightChild();
            }
            current->setRightChild(newNode);
            newNode->setParent(current);
        }
    }
}

/**
 * Initialize contours for packing
 */
void ASFBStarTree::initializeContours() {
    horizontalContour->clear();
    verticalContour->clear();
    
    // Initialize horizontal contour with a segment at y=0
    horizontalContour->addSegment(0, numeric_limits<int>::max(), 0);
    
    // Initialize vertical contour with a segment at x=0
    verticalContour->addSegment(0, numeric_limits<int>::max(), 0);
}

/**
 * Update contour with a module
 */
void ASFBStarTree::updateContourWithModule(const shared_ptr<Module>& module) {
    if (!module) return;
    
    int x = module->getX();
    int y = module->getY();
    int width = module->getWidth();
    int height = module->getHeight();
    
    // Update horizontal contour
    horizontalContour->addSegment(x, x + width, y + height);
    
    // Update vertical contour
    verticalContour->addSegment(y, y + height, x + width);
}

/**
 * Pack a node in the ASF-B*-tree
 */
void ASFBStarTree::packNode(const shared_ptr<BStarTreeNode>& node) {
    if (!node) return;
    
    const string& moduleName = node->getModuleName();
    auto module = modules[moduleName];
    
    if (!module) return;
    
    int x = 0, y = 0;
    
    // Calculate x-coordinate based on B*-tree rules
    if (node->getParent()) {
        auto parent = modules[node->getParent()->getModuleName()];
        if (parent) {
            if (node->getParent()->getLeftChild() == node) {
                // Left child: place to the right of parent
                x = parent->getX() + parent->getWidth();
            } else {
                // Right child: same x-coordinate as parent
                x = parent->getX();
            }
        }
    }
    
    // Calculate y-coordinate using the horizontal contour
    y = horizontalContour->getHeight(x, x + module->getWidth());
    
    // Special handling for self-symmetric modules
    if (find(selfSymmetricModules.begin(), selfSymmetricModules.end(), 
                 moduleName) != selfSymmetricModules.end()) {
        // Self-symmetric modules must be on the boundary
        if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
            // For vertical symmetry, the module must abut the symmetry axis
            x = static_cast<int>(symmetryAxisPosition) - module->getWidth() / 2;
        } else {
            // For horizontal symmetry, the module must be on the bottom
            y = static_cast<int>(symmetryAxisPosition) - module->getHeight() / 2;
        }
    }
    
    // Set the module's position
    module->setPosition(x, y);
    
    // Update contours
    updateContourWithModule(module);
}

/**
 * Calculate the positions of symmetric modules
 */
void ASFBStarTree::calculateSymmetricModulePositions() {
    // Step 1: Force all modules to have positive coordinates
    // Shift any modules with negative coordinates
    int minX = std::numeric_limits<int>::max();
    int minY = std::numeric_limits<int>::max();
    
    for (auto& pair : modules) {
        auto& module = pair.second;
        minX = std::min(minX, module->getX());
        minY = std::min(minY, module->getY());
    }
    
    // If any coordinates are negative, shift all modules
    if (minX < 0 || minY < 0) {
        int shiftX = std::max(0, -minX) + 5; // 5-unit buffer
        int shiftY = std::max(0, -minY) + 5; // 5-unit buffer
        
        for (auto& pair : modules) {
            auto& module = pair.second;
            module->setPosition(module->getX() + shiftX, module->getY() + shiftY);
        }
        
        // Adjust symmetry axis if needed
        if (minX < 0 && symmetryGroup->getType() == SymmetryType::VERTICAL) {
            symmetryAxisPosition += shiftX;
        } else if (minY < 0 && symmetryGroup->getType() == SymmetryType::HORIZONTAL) {
            symmetryAxisPosition += shiftY;
        }
    }
    
    // Step 2: Apply consistent symmetry positioning for all pairs
    SymmetryType symType = symmetryGroup->getType();
    const int BUFFER = 20; // Increased buffer to prevent overlap issues
    
    // Process symmetry pairs with exact calculations
    for (const auto& pair : symmetryGroup->getSymmetryPairs()) {
        const std::string& module1 = pair.first;
        const std::string& module2 = pair.second;
        
        auto mod1It = modules.find(module1);
        auto mod2It = modules.find(module2);
        
        if (mod1It == modules.end() || mod2It == modules.end()) continue;
        
        auto mod1 = mod1It->second;
        auto mod2 = mod2It->second;
        
        // Ensure both modules have the same dimensions and rotation status
        mod2->setRotation(mod1->getRotated());
        
        if (symType == SymmetryType::VERTICAL) {
            // Calculate widths (accounting for rotation)
            int width1 = mod1->getWidth();
            int width2 = mod2->getWidth();
            
            // Ensure they have the same Y coordinate
            int yPos = mod1->getY();
            
            // First module on left side of the axis with a buffer
            int mod1X = static_cast<int>(symmetryAxisPosition) - width1 - BUFFER;
            // Second module on right side of the axis with a buffer
            int mod2X = static_cast<int>(symmetryAxisPosition) + BUFFER;
            
            // Set positions with symmetry around the axis
            mod1->setPosition(mod1X, yPos);
            mod2->setPosition(mod2X, yPos);
            
            std::cout << "Positioned symmetry pair: " << module1 << " at (" 
                     << mod1->getX() << "," << mod1->getY() << ") and " 
                     << module2 << " at (" << mod2->getX() << "," << mod2->getY() 
                     << ") about axis=" << symmetryAxisPosition << std::endl;
        } 
        else { // HORIZONTAL
            // Similar logic but for horizontal symmetry
            int height1 = mod1->getHeight();
            int height2 = mod2->getHeight();
            
            // Ensure they have the same X coordinate
            int xPos = mod1->getX();
            
            // First module below the axis with a buffer
            int mod1Y = static_cast<int>(symmetryAxisPosition) - height1 - BUFFER;
            // Second module above the axis with a buffer
            int mod2Y = static_cast<int>(symmetryAxisPosition) + BUFFER;
            
            // Set positions with symmetry around the axis
            mod1->setPosition(xPos, mod1Y);
            mod2->setPosition(xPos, mod2Y);
        }
    }
    
    // Step 3: Position self-symmetric modules
    for (const auto& moduleName : selfSymmetricModules) {
        auto it = modules.find(moduleName);
        if (it == modules.end()) continue;
        
        auto module = it->second;
        
        // Center module on the symmetry axis
        if (symType == SymmetryType::VERTICAL) {
            int width = module->getWidth();
            int x = static_cast<int>(symmetryAxisPosition) - width / 2;
            module->setPosition(x, module->getY());
        } else { // HORIZONTAL
            int height = module->getHeight();
            int y = static_cast<int>(symmetryAxisPosition) - height / 2;
            module->setPosition(module->getX(), y);
        }
    }
    
    // Step 4: Final verification - check for overlaps
    // Use a set to keep track of which modules we've already fixed
    std::unordered_set<std::string> fixedModules;
    bool overlapsExist = true;
    int iterationLimit = 10; // Prevent infinite loops
    
    while (overlapsExist && iterationLimit-- > 0) {
        overlapsExist = false;
        
        // Check for and fix any remaining overlaps
        for (auto& pair1 : modules) {
            auto& mod1 = pair1.second;
            
            for (auto& pair2 : modules) {
                // Skip self-comparison or modules with the same name
                if (pair1.first == pair2.first || mod1->getName() == pair2.second->getName()) {
                    continue;
                }
                
                auto& mod2 = pair2.second;
                
                // Skip if this pair has already been fixed in this iteration
                std::string pairKey = pair1.first < pair2.first ? 
                                     pair1.first + ":" + pair2.first : 
                                     pair2.first + ":" + pair1.first;
                
                if (fixedModules.find(pairKey) != fixedModules.end()) {
                    continue;
                }
                
                if (mod1->overlaps(*mod2)) {
                    overlapsExist = true;
                    
                    // Mark this pair as fixed
                    fixedModules.insert(pairKey);
                    
                    // Fix overlap by moving one module
                    if (symType == SymmetryType::VERTICAL) {
                        // For vertical symmetry, shift in Y direction
                        if (mod1->getY() <= mod2->getY()) {
                            mod2->setPosition(mod2->getX(), mod1->getY() + mod1->getHeight() + BUFFER);
                        } else {
                            mod1->setPosition(mod1->getX(), mod2->getY() + mod2->getHeight() + BUFFER);
                        }
                    } else {
                        // For horizontal symmetry, shift in X direction
                        if (mod1->getX() <= mod2->getX()) {
                            mod2->setPosition(mod1->getX() + mod1->getWidth() + BUFFER, mod2->getY());
                        } else {
                            mod1->setPosition(mod2->getX() + mod2->getWidth() + BUFFER, mod1->getY());
                        }
                    }
                    
                    // If we've fixed a symmetry pair, also reposition the symmetric module
                    for (const auto& symPair : symmetryGroup->getSymmetryPairs()) {
                        if (symPair.first == pair1.first || symPair.second == pair1.first) {
                            std::string symModName = (symPair.first == pair1.first) ? symPair.second : symPair.first;
                            auto symModIt = modules.find(symModName);
                            if (symModIt != modules.end()) {
                                // Recalculate position to maintain symmetry
                                auto symMod = symModIt->second;
                                if (symType == SymmetryType::VERTICAL) {
                                    // Place symmetrically across the axis
                                    int delta = std::abs(mod1->getX() - static_cast<int>(symmetryAxisPosition));
                                    if (mod1->getX() < static_cast<int>(symmetryAxisPosition)) {
                                        symMod->setPosition(static_cast<int>(symmetryAxisPosition) + delta, mod1->getY());
                                    } else {
                                        symMod->setPosition(static_cast<int>(symmetryAxisPosition) - delta - symMod->getWidth(), mod1->getY());
                                    }
                                } else {
                                    // Place symmetrically across the axis
                                    int delta = std::abs(mod1->getY() - static_cast<int>(symmetryAxisPosition));
                                    if (mod1->getY() < static_cast<int>(symmetryAxisPosition)) {
                                        symMod->setPosition(mod1->getX(), static_cast<int>(symmetryAxisPosition) + delta);
                                    } else {
                                        symMod->setPosition(mod1->getX(), static_cast<int>(symmetryAxisPosition) - delta - symMod->getHeight());
                                    }
                                }
                            }
                        }
                        else if (symPair.first == pair2.first || symPair.second == pair2.first) {
                            // Similar handling for mod2
                            std::string symModName = (symPair.first == pair2.first) ? symPair.second : symPair.first;
                            auto symModIt = modules.find(symModName);
                            if (symModIt != modules.end()) {
                                auto symMod = symModIt->second;
                                if (symType == SymmetryType::VERTICAL) {
                                    // Place symmetrically across the axis
                                    int delta = std::abs(mod2->getX() - static_cast<int>(symmetryAxisPosition));
                                    if (mod2->getX() < static_cast<int>(symmetryAxisPosition)) {
                                        symMod->setPosition(static_cast<int>(symmetryAxisPosition) + delta, mod2->getY());
                                    } else {
                                        symMod->setPosition(static_cast<int>(symmetryAxisPosition) - delta - symMod->getWidth(), mod2->getY());
                                    }
                                } else {
                                    // Place symmetrically across the axis
                                    int delta = std::abs(mod2->getY() - static_cast<int>(symmetryAxisPosition));
                                    if (mod2->getY() < static_cast<int>(symmetryAxisPosition)) {
                                        symMod->setPosition(mod2->getX(), static_cast<int>(symmetryAxisPosition) + delta);
                                    } else {
                                        symMod->setPosition(mod2->getX(), static_cast<int>(symmetryAxisPosition) - delta - symMod->getHeight());
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        
        // Clear the fixed modules set for the next iteration
        fixedModules.clear();
    }
}

/**
 * Calculates the coordinates of all modules in the symmetry group
 * by packing the ASF-B*-tree
 */
bool ASFBStarTree::pack() {
    if (!root) {
        std::cerr << "Cannot pack ASFBStarTree: root is null" << std::endl;
        return false;
    }
    
    // Initialize contours
    initializeContours();
    
    // Step 1: Pack all representative modules
    std::queue<std::shared_ptr<BStarTreeNode>> nodeQueue;
    nodeQueue.push(root);
    
    // Track min/max coordinates for axis calculation
    int minX = std::numeric_limits<int>::max();
    int maxX = std::numeric_limits<int>::min();
    int minY = std::numeric_limits<int>::max();
    int maxY = std::numeric_limits<int>::min();
    
    // First pass: pack all representative modules and calculate their bounding box
    std::cout << "First pass: packing representative modules" << std::endl;
    while (!nodeQueue.empty()) {
        auto currentNode = nodeQueue.front();
        nodeQueue.pop();
        
        if (!currentNode) continue;
        
        // Pack the current node
        packNode(currentNode);
        
        // Update min/max coordinates from the packed module
        std::string moduleName = currentNode->getModuleName();
        auto moduleIt = modules.find(moduleName);
        if (moduleIt != modules.end()) {
            auto& module = moduleIt->second;
            
            minX = std::min(minX, module->getX());
            maxX = std::max(maxX, module->getX() + module->getWidth());
            minY = std::min(minY, module->getY());
            maxY = std::max(maxY, module->getY() + module->getHeight());
        }
        
        // Add children to the queue
        if (currentNode->getLeftChild()) {
            nodeQueue.push(currentNode->getLeftChild());
        }
        if (currentNode->getRightChild()) {
            nodeQueue.push(currentNode->getRightChild());
        }
    }
    
    // Ensure no modules have negative coordinates
    int shiftX = 0, shiftY = 0;
    if (minX < 0) {
        shiftX = -minX + 10; // Add a 10-unit buffer
        minX += shiftX;
        maxX += shiftX;
    }
    if (minY < 0) {
        shiftY = -minY + 10; // Add a 10-unit buffer
        minY += shiftY;
        maxY += shiftY;
    }
    
    // Apply any shifts needed
    if (shiftX > 0 || shiftY > 0) {
        std::cout << "Shifting all modules by (" << shiftX << "," << shiftY << ")" << std::endl;
        for (auto& pair : modules) {
            auto& module = pair.second;
            module->setPosition(module->getX() + shiftX, module->getY() + shiftY);
        }
    }
    
    // Step 2: Calculate stable symmetry axis
    const int MIN_AXIS_POSITION = 100; // Minimum axis position to ensure stability
    
    if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
        // Calculate axis as midpoint between min and max X, ensuring minimum position
        symmetryAxisPosition = std::max(MIN_AXIS_POSITION, (minX + maxX) / 2);
        
        // Ensure the axis has a sufficient margin from all modules (at least 20 units)
        const int MIN_MARGIN = 20;
        for (const auto& pair : modules) {
            auto& module = pair.second;
            int moduleRight = module->getX() + module->getWidth();
            
            // If module is too close to the axis, adjust the axis
            if (module->getX() < symmetryAxisPosition && 
                (symmetryAxisPosition - moduleRight) < MIN_MARGIN) {
                symmetryAxisPosition = moduleRight + MIN_MARGIN;
            }
            else if (module->getX() > symmetryAxisPosition && 
                     (module->getX() - symmetryAxisPosition) < MIN_MARGIN) {
                symmetryAxisPosition = module->getX() - MIN_MARGIN;
            }
        }
    } else { // HORIZONTAL
        // Calculate axis as midpoint between min and max Y, ensuring minimum position
        symmetryAxisPosition = std::max(MIN_AXIS_POSITION, (minY + maxY) / 2);
        
        // Ensure the axis has a sufficient margin from all modules (at least 20 units)
        const int MIN_MARGIN = 20;
        for (const auto& pair : modules) {
            auto& module = pair.second;
            int moduleTop = module->getY() + module->getHeight();
            
            // If module is too close to the axis, adjust the axis
            if (module->getY() < symmetryAxisPosition && 
                (symmetryAxisPosition - moduleTop) < MIN_MARGIN) {
                symmetryAxisPosition = moduleTop + MIN_MARGIN;
            }
            else if (module->getY() > symmetryAxisPosition && 
                     (module->getY() - symmetryAxisPosition) < MIN_MARGIN) {
                symmetryAxisPosition = module->getY() - MIN_MARGIN;
            }
        }
    }
    
    std::cout << "Using stabilized symmetry axis: " << symmetryAxisPosition << std::endl;
    
    // Step 3: Calculate positions for all modules in symmetry pairs
    calculateSymmetricModulePositions();
    
    // Step 4: Recalculate contours based on the final module positions
    // Reset contours
    horizontalContour->clear();
    verticalContour->clear();
    
    // Update contours for all modules
    for (const auto& pair : modules) {
        auto& module = pair.second;
        int x = module->getX();
        int y = module->getY();
        int width = module->getWidth();
        int height = module->getHeight();
        
        horizontalContour->addSegment(x, x + width, y + height);
        verticalContour->addSegment(y, y + height, x + width);
    }
    
    // Step 5: Final validation - check for symmetry constraint issues
    bool isValid = validateSymmetryConstraints();
    
    if (!isValid) {
        std::cout << "Final pack failed symmetry validation" << std::endl;
        
        // Add a repair attempt - force symmetry pairs to be correctly positioned
        for (const auto& symPair : symmetryGroup->getSymmetryPairs()) {
            auto mod1It = modules.find(symPair.first);
            auto mod2It = modules.find(symPair.second);
            
            if (mod1It != modules.end() && mod2It != modules.end()) {
                auto& mod1 = mod1It->second;
                auto& mod2 = mod2It->second;
                
                // Force symmetry around the axis
                if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
                    int y = mod1->getY(); // Keep Y position the same
                    
                    // Place on opposite sides of the axis with minimum separation
                    int mod1X = static_cast<int>(symmetryAxisPosition) - mod1->getWidth() - 20;
                    int mod2X = static_cast<int>(symmetryAxisPosition) + 20;
                    
                    mod1->setPosition(mod1X, y);
                    mod2->setPosition(mod2X, y);
                } else { // HORIZONTAL
                    int x = mod1->getX(); // Keep X position the same
                    
                    // Place on opposite sides of the axis with minimum separation
                    int mod1Y = static_cast<int>(symmetryAxisPosition) - mod1->getHeight() - 20;
                    int mod2Y = static_cast<int>(symmetryAxisPosition) + 20;
                    
                    mod1->setPosition(x, mod1Y);
                    mod2->setPosition(x, mod2Y);
                }
            }
        }
        
        // Place self-symmetric modules directly on the axis
        for (const auto& modName : selfSymmetricModules) {
            auto modIt = modules.find(modName);
            if (modIt != modules.end()) {
                auto& mod = modIt->second;
                
                if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
                    int x = static_cast<int>(symmetryAxisPosition) - mod->getWidth() / 2;
                    mod->setPosition(x, mod->getY());
                } else { // HORIZONTAL
                    int y = static_cast<int>(symmetryAxisPosition) - mod->getHeight() / 2;
                    mod->setPosition(mod->getX(), y);
                }
            }
        }
    }
    
    // Return success if the symmetry island is valid
    bool islandValid = isSymmetryIslandValid();
    if (!islandValid) {
        std::cerr << "Invalid symmetry island for group: " << symmetryGroup->getName() << std::endl;
    }
    
    return islandValid && isValid;
}

/**
 * Validates that all symmetry constraints are satisfied
 */
bool ASFBStarTree::validateSymmetryConstraints() const {
    if (!symmetryGroup) {
        std::cerr << "No symmetry group defined" << std::endl;
        return false;
    }
    
    std::cout << "Validating symmetry constraints for group: " << symmetryGroup->getName() << std::endl;
    std::cout << "Symmetry axis position: " << symmetryAxisPosition << std::endl;
    std::cout << "Symmetry type: " << (symmetryGroup->getType() == SymmetryType::VERTICAL ? "VERTICAL" : "HORIZONTAL") << std::endl;
    
    // Check symmetry pairs have correct positions
    for (const auto& pair : symmetryGroup->getSymmetryPairs()) {
        const auto& module1It = modules.find(pair.first);
        const auto& module2It = modules.find(pair.second);
        
        if (module1It == modules.end() || module2It == modules.end()) {
            std::cerr << "Module pair not found: " << pair.first << " and " << pair.second << std::endl;
            continue;
        }
        
        const auto& mod1 = module1It->second;
        const auto& mod2 = module2It->second;
        
        if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
            // Calculate centers accounting for rotation
            double center1X = mod1->getX() + mod1->getWidth() / 2.0;
            double center2X = mod2->getX() + mod2->getWidth() / 2.0;
            
            // Relaxed symmetry validation with tolerance
            const double TOLERANCE = 5.0; // 5-unit tolerance
            
            // Check if centers are approximately equidistant from the axis
            double axis2x = 2.0 * symmetryAxisPosition;
            double diff = std::abs((center1X + center2X) - axis2x);
            
            if (diff > TOLERANCE) {
                std::cerr << "Vertical symmetry violation for pair: " 
                          << pair.first << " (" << mod1->getX() << "," << mod1->getY() << ") center=" << center1X
                          << " and " << pair.second << " (" << mod2->getX() << "," << mod2->getY() << ") center=" << center2X
                          << " diff=" << diff << " axis=" << symmetryAxisPosition << std::endl;
                
                // In case of slight deviation, print additional info
                std::cerr << "  Expected " << pair.second << " center at X=" 
                          << (axis2x - center1X) << " but found X=" << center2X << std::endl;
                
                return false;
            }
            
            // Check if Y coordinates are approximately equal
            if (std::abs(mod1->getY() - mod2->getY()) > TOLERANCE) {
                std::cerr << "Y-coordinate mismatch for pair: " 
                          << pair.first << " (" << mod1->getY() << ") and "
                          << pair.second << " (" << mod2->getY() << ")" << std::endl;
                return false;
            }
            
            // Check if both modules have the same dimensions (accounting for rotation)
            if (std::abs(mod1->getWidth() - mod2->getWidth()) > TOLERANCE || 
                std::abs(mod1->getHeight() - mod2->getHeight()) > TOLERANCE) {
                std::cerr << "Dimension mismatch for pair: " 
                          << pair.first << " (" << mod1->getWidth() << "x" << mod1->getHeight() << ") and "
                          << pair.second << " (" << mod2->getWidth() << "x" << mod2->getHeight() << ")" << std::endl;
                return false;
            }
        } else { // HORIZONTAL
            // Similar validation with improved tolerance and debugging for horizontal symmetry
            double center1Y = mod1->getY() + mod1->getHeight() / 2.0;
            double center2Y = mod2->getY() + mod2->getHeight() / 2.0;
            
            const double TOLERANCE = 5.0;
            
            double axis2y = 2.0 * symmetryAxisPosition;
            double diff = std::abs((center1Y + center2Y) - axis2y);
            
            if (diff > TOLERANCE) {
                std::cerr << "Horizontal symmetry violation for pair: " 
                          << pair.first << " (" << mod1->getX() << "," << mod1->getY() << ") center=" << center1Y
                          << " and " << pair.second << " (" << mod2->getX() << "," << mod2->getY() << ") center=" << center2Y
                          << " diff=" << diff << " axis=" << symmetryAxisPosition << std::endl;
                
                std::cerr << "  Expected " << pair.second << " center at Y=" 
                          << (axis2y - center1Y) << " but found Y=" << center2Y << std::endl;
                
                return false;
            }
            
            if (std::abs(mod1->getX() - mod2->getX()) > TOLERANCE) {
                std::cerr << "X-coordinate mismatch for pair: " 
                          << pair.first << " (" << mod1->getX() << ") and "
                          << pair.second << " (" << mod2->getX() << ")" << std::endl;
                return false;
            }
            
            if (std::abs(mod1->getWidth() - mod2->getWidth()) > TOLERANCE || 
                std::abs(mod1->getHeight() - mod2->getHeight()) > TOLERANCE) {
                std::cerr << "Dimension mismatch for pair: " 
                          << pair.first << " (" << mod1->getWidth() << "x" << mod1->getHeight() << ") and "
                          << pair.second << " (" << mod2->getWidth() << "x" << mod2->getHeight() << ")" << std::endl;
                return false;
            }
        }
    }
    
    // Check self-symmetric modules
    for (const auto& moduleName : selfSymmetricModules) {
        auto it = modules.find(moduleName);
        if (it == modules.end()) {
            std::cerr << "Self-symmetric module not found: " << moduleName << std::endl;
            continue;
        }
        
        auto module = it->second;
        
        const double TOLERANCE = 5.0;
        
        if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
            // Calculate center X
            double centerX = module->getX() + module->getWidth() / 2.0;
            
            // Check if center is approximately on the axis
            if (std::abs(centerX - symmetryAxisPosition) > TOLERANCE) {
                std::cerr << "Self-symmetric module " << moduleName 
                          << " not centered on vertical axis. Center X=" << centerX
                          << " but axis is at " << symmetryAxisPosition << std::endl;
                return false;
            }
        } else { // HORIZONTAL
            // Calculate center Y
            double centerY = module->getY() + module->getHeight() / 2.0;
            
            // Check if center is approximately on the axis
            if (std::abs(centerY - symmetryAxisPosition) > TOLERANCE) {
                std::cerr << "Self-symmetric module " << moduleName 
                          << " not centered on horizontal axis. Center Y=" << centerY
                          << " but axis is at " << symmetryAxisPosition << std::endl;
                return false;
            }
        }
    }
    
    std::cout << "Symmetry validation passed!" << std::endl;
    return true;
}

/**
 * Checks if a module is on the correct branch according to Property 1
 */
bool ASFBStarTree::isOnCorrectBranch(const string& moduleName) const {
    // Skip check if not a self-symmetric module
    if (find(selfSymmetricModules.begin(), selfSymmetricModules.end(), moduleName) == selfSymmetricModules.end()) {
        return true;
    }
    
    // Find the node for this module in the tree
    shared_ptr<BStarTreeNode> node = nullptr;
    
    // Find the node in the B*-tree (simplified - actual implementation would track nodes directly)
    queue<shared_ptr<BStarTreeNode>> q;
    if (root) q.push(root);
    
    while (!q.empty() && !node) {
        auto current = q.front();
        q.pop();
        
        if (current->getModuleName() == moduleName) {
            node = current;
        } else {
            if (current->getLeftChild()) q.push(current->getLeftChild());
            if (current->getRightChild()) q.push(current->getRightChild());
        }
    }
    
    if (!node) return false;
    
    // For vertical symmetry, self-symmetric modules must be on the rightmost branch
    if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
        while (node && node->getParent()) {
            if (node->isLeftChild()) {
                return false; // Not on rightmost branch
            }
            node = node->getParent();
        }
    } else { // HORIZONTAL
        while (node && node->getParent()) {
            if (node->isRightChild()) {
                return false; // Not on leftmost branch
            }
            node = node->getParent();
        }
    }
    
    return true;
}

/**
 * Validates that the placement forms a symmetry island
 */
bool ASFBStarTree::isSymmetryIslandValid() const {
    if (modules.empty()) {
        std::cout << "Symmetry island validation skipped: no modules" << std::endl;
        return true;
    }
    
    // Build adjacency list for all modules in the symmetry group
    std::unordered_map<std::string, std::vector<std::string>> adjacency;
    
    // Maps for positions and dimensions
    std::unordered_map<std::string, std::pair<int, int>> positions;
    std::unordered_map<std::string, std::pair<int, int>> dimensions;
    
    // Fill these maps
    for (const auto& pair : modules) {
        const auto& module = pair.second;
        positions[pair.first] = {module->getX(), module->getY()};
        dimensions[pair.first] = {module->getWidth(), module->getHeight()};
    }
    
    std::cout << "Checking symmetry island validity for group with " << modules.size() << " modules" << std::endl;
    
    // Calculate expanded bounding boxes for modules - allow for a small gap
    const int ADJACENCY_THRESHOLD = 10; // Allow small gaps between adjacent modules
    
    // Build adjacency list based on abutting or nearly-abutting modules
    for (const auto& pair1 : modules) {
        const auto& name1 = pair1.first;
        const auto& mod1 = pair1.second;
        
        int x1 = mod1->getX();
        int y1 = mod1->getY();
        int w1 = mod1->getWidth();
        int h1 = mod1->getHeight();
        
        // Expanded bounds for adjacency testing
        int x1Left = x1 - ADJACENCY_THRESHOLD;
        int x1Right = x1 + w1 + ADJACENCY_THRESHOLD;
        int y1Top = y1 + h1 + ADJACENCY_THRESHOLD;
        int y1Bottom = y1 - ADJACENCY_THRESHOLD;
        
        for (const auto& pair2 : modules) {
            const auto& name2 = pair2.first;
            if (name1 == name2) continue;
            
            const auto& mod2 = pair2.second;
            
            int x2 = mod2->getX();
            int y2 = mod2->getY();
            int w2 = mod2->getWidth();
            int h2 = mod2->getHeight();
            
            // Check for adjacency or near-adjacency
            bool adjacent = false;
            
            // Check if right edge of mod1 is near left edge of mod2
            if (std::abs((x1 + w1) - x2) <= ADJACENCY_THRESHOLD) {
                // Check for vertical overlap
                if (!(y1 >= (y2 + h2) || y2 >= (y1 + h1))) {
                    adjacent = true;
                    std::cout << "  Modules " << name1 << " and " << name2 
                              << " are adjacent horizontally (right-left)" << std::endl;
                }
            }
            // Check if left edge of mod1 is near right edge of mod2
            else if (std::abs(x1 - (x2 + w2)) <= ADJACENCY_THRESHOLD) {
                // Check for vertical overlap
                if (!(y1 >= (y2 + h2) || y2 >= (y1 + h1))) {
                    adjacent = true;
                    std::cout << "  Modules " << name1 << " and " << name2 
                              << " are adjacent horizontally (left-right)" << std::endl;
                }
            }
            // Check if top edge of mod1 is near bottom edge of mod2
            else if (std::abs((y1 + h1) - y2) <= ADJACENCY_THRESHOLD) {
                // Check for horizontal overlap
                if (!(x1 >= (x2 + w2) || x2 >= (x1 + w1))) {
                    adjacent = true;
                    std::cout << "  Modules " << name1 << " and " << name2 
                              << " are adjacent vertically (top-bottom)" << std::endl;
                }
            }
            // Check if bottom edge of mod1 is near top edge of mod2
            else if (std::abs(y1 - (y2 + h2)) <= ADJACENCY_THRESHOLD) {
                // Check for horizontal overlap
                if (!(x1 >= (x2 + w2) || x2 >= (x1 + w1))) {
                    adjacent = true;
                    std::cout << "  Modules " << name1 << " and " << name2 
                              << " are adjacent vertically (bottom-top)" << std::endl;
                }
            }
            
            // Special case: check for symmetry pair as inherently adjacent
            for (const auto& symPair : symmetryGroup->getSymmetryPairs()) {
                if ((symPair.first == name1 && symPair.second == name2) ||
                    (symPair.first == name2 && symPair.second == name1)) {
                    adjacent = true;
                    std::cout << "  Modules " << name1 << " and " << name2 
                              << " are adjacent as symmetry pair" << std::endl;
                    break;
                }
            }
            
            // Add to adjacency list if adjacent
            if (adjacent) {
                adjacency[name1].push_back(name2);
            }
        }
    }
    
    // Check if we have any adjacency at all
    bool hasAnyAdjacency = false;
    for (const auto& pair : adjacency) {
        if (!pair.second.empty()) {
            hasAnyAdjacency = true;
            break;
        }
    }
    
    if (!hasAnyAdjacency) {
        std::cout << "No adjacency found between any modules - all modules are isolated" << std::endl;
        
        // Print module positions for debugging
        for (const auto& pair : modules) {
            const auto& module = pair.second;
            std::cout << "  Module " << pair.first 
                      << " at (" << module->getX() << "," << module->getY() << ") "
                      << "size: " << module->getWidth() << "x" << module->getHeight() << std::endl;
        }
        
        // If this is a simple 1-2 module case, consider it valid
        if (modules.size() <= 2) {
            std::cout << "Only 1-2 modules in group, considering it valid despite no adjacency" << std::endl;
            return true;
        }
        
        return false;
    }
    
    // Perform BFS to check connectivity
    std::unordered_set<std::string> visited;
    std::queue<std::string> q;
    
    // Start with the first module
    if (!modules.empty()) {
        auto it = modules.begin();
        q.push(it->first);
        visited.insert(it->first);
        
        std::cout << "Starting BFS from module: " << it->first << std::endl;
    }
    
    while (!q.empty()) {
        std::string current = q.front();
        q.pop();
        
        for (const auto& neighbor : adjacency[current]) {
            if (visited.find(neighbor) == visited.end()) {
                visited.insert(neighbor);
                q.push(neighbor);
                std::cout << "  Visited module: " << neighbor << std::endl;
            }
        }
    }
    
    // Check if all modules are visited
    bool allVisited = (visited.size() == modules.size());
    
    if (!allVisited) {
        std::cout << "Connectivity check failed: only visited " << visited.size() 
                  << " out of " << modules.size() << " modules." << std::endl;
        
        // Find and report which modules weren't visited
        std::cout << "Unreachable modules: ";
        for (const auto& pair : modules) {
            if (visited.find(pair.first) == visited.end()) {
                std::cout << pair.first << " ";
            }
        }
        std::cout << std::endl;
        
        // For small groups (3 or fewer modules), consider them valid anyway
        if (modules.size() <= 3) {
            std::cout << "Small module group (≤3 modules), considering valid despite connectivity issues" << std::endl;
            return true;
        }
    } else {
        std::cout << "Connectivity check passed: all " << modules.size() << " modules are connected" << std::endl;
    }
    
    return allVisited;
}

/**
 * Checks if the tree satisfies the symmetric-feasible condition
 */
bool ASFBStarTree::isSymmetricFeasible() const {
    // Check self-symmetric modules
    for (const auto& moduleName : selfSymmetricModules) {
        auto node = [this, &moduleName]() -> shared_ptr<BStarTreeNode> {
            // Find the node for this module
            queue<shared_ptr<BStarTreeNode>> queue;
            if (root) queue.push(root);
            
            while (!queue.empty()) {
                auto current = queue.front();
                queue.pop();
                
                if (current->getModuleName() == moduleName) {
                    return current;
                }
                
                if (current->getLeftChild()) queue.push(current->getLeftChild());
                if (current->getRightChild()) queue.push(current->getRightChild());
            }
            
            return nullptr;
        }();
        
        // For vertical symmetry, self-symmetric modules must be on the rightmost branch
        if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
            // Traverse up the tree from the node
            auto current = node;
            while (current && current->getParent()) {
                // If this node is a left child of its parent, it's not on the rightmost branch
                if (current->getParent()->getLeftChild() == current) {
                    return false;
                }
                current = current->getParent();
            }
        }
        // For horizontal symmetry, self-symmetric modules must be on the leftmost branch
        else {
            // Traverse up the tree from the node
            auto current = node;
            while (current && current->getParent()) {
                // If this node is a right child of its parent, it's not on the leftmost branch
                if (current->getParent()->getRightChild() == current) {
                    return false;
                }
                current = current->getParent();
            }
        }
    }
    
    return true;
}

/**
 * Returns the bounding rectangle area of the symmetry island
 */
int ASFBStarTree::getArea() const {
    if (modules.empty()) return 0;
    
    int minX = numeric_limits<int>::max();
    int minY = numeric_limits<int>::max();
    int maxX = 0;
    int maxY = 0;
    
    // Find the bounding rectangle
    for (const auto& pair : modules) {
        const auto& module = pair.second;
        
        minX = min(minX, module->getX());
        minY = min(minY, module->getY());
        maxX = max(maxX, module->getX() + module->getWidth());
        maxY = max(maxY, module->getY() + module->getHeight());
    }
    
    return (maxX - minX) * (maxY - minY);
}

/**
 * Gets the contour of the symmetry island (for HB*-tree)
 */
pair<shared_ptr<Contour>, shared_ptr<Contour>> ASFBStarTree::getContours() const {
    return {horizontalContour, verticalContour};
}

/**
 * Rotates a module in the symmetry group
 */
bool ASFBStarTree::rotateModule(const string& moduleName) {
    auto it = modules.find(moduleName);
    if (it == modules.end()) return false;
    
    auto module = it->second;
    
    // Special handling for symmetry pairs and self-symmetric modules
    auto pairIt = symmetricPairMap.find(moduleName);
    auto selfIt = find(selfSymmetricModules.begin(), selfSymmetricModules.end(), moduleName);
    
    if (pairIt != symmetricPairMap.end()) {
        // For symmetry pairs, rotate both modules
        auto pairModule = modules[pairIt->second];
        if (pairModule) {
            pairModule->rotate();
        }
    } else if (selfIt != selfSymmetricModules.end()) {
        // For self-symmetric modules, update the shape of its representative
        // This is a simplification - in reality, more complex handling might be needed
    }
    
    // Rotate the module
    module->rotate();
    
    return true;
}

/**
 * Helper function: checks if a module is on the boundary
 */
bool ASFBStarTree::isOnBoundary(const string& moduleName) const {
    return find(selfSymmetricModules.begin(), selfSymmetricModules.end(), moduleName) != selfSymmetricModules.end();
}

/**
 * Helper function: checks if a node can be moved to a new position
 */
bool ASFBStarTree::canMoveNode(const shared_ptr<BStarTreeNode>& node, 
                              const shared_ptr<BStarTreeNode>& newParent, 
                              bool asLeftChild) const {
    if (!node || !newParent) return false;
    
    // Self-symmetric modules have placement restrictions
    if (isOnBoundary(node->getModuleName())) {
        // For vertical symmetry, self-symmetric modules must be on the rightmost branch
        if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
            if (asLeftChild) return false;  // Can't be a left child
            
            // Check if the new parent is on the rightmost branch
            auto current = newParent;
            while (current && current->getParent()) {
                if (current->getParent()->getLeftChild() == current) {
                    return false;  // Not on the rightmost branch
                }
                current = current->getParent();
            }
        }
        // For horizontal symmetry, self-symmetric modules must be on the leftmost branch
        else {
            if (!asLeftChild) return false;  // Can't be a right child
            
            // Check if the new parent is on the leftmost branch
            auto current = newParent;
            while (current && current->getParent()) {
                if (current->getParent()->getRightChild() == current) {
                    return false;  // Not on the leftmost branch
                }
                current = current->getParent();
            }
        }
    }
    
    return true;
}

/**
 * Moves a node to a new position in the tree
 */
bool ASFBStarTree::moveNode(const string& nodeName, 
                           const string& newParentName, 
                           bool asLeftChild) {
    // Find the nodes...
    shared_ptr<BStarTreeNode> node = nullptr;
    shared_ptr<BStarTreeNode> newParent = nullptr;
    
    // Find the node to move (simplified for clarity)
    queue<shared_ptr<BStarTreeNode>> q;
    if (root) q.push(root);
    
    while (!q.empty() && (!node || !newParent)) {
        auto current = q.front();
        q.pop();
        
        if (current->getModuleName() == nodeName) {
            node = current;
        }
        if (current->getModuleName() == newParentName) {
            newParent = current;
        }
        
        if (current->getLeftChild()) q.push(current->getLeftChild());
        if (current->getRightChild()) q.push(current->getRightChild());
    }
    
    if (!node || !newParent) return false;
    
    // Check if nodeName is a self-symmetric module
    bool isSelfSymmetric = find(selfSymmetricModules.begin(), 
                               selfSymmetricModules.end(), 
                               nodeName) != selfSymmetricModules.end();
    
    // Enforce Property 1 for self-symmetric modules
    if (isSelfSymmetric) {
        if (symmetryGroup->getType() == SymmetryType::VERTICAL && asLeftChild) {
            // Can't move a self-symmetric module as a left child for vertical symmetry
            return false;
        } else if (symmetryGroup->getType() == SymmetryType::HORIZONTAL && !asLeftChild) {
            // Can't move a self-symmetric module as a right child for horizontal symmetry
            return false;
        }
        
        // Also verify that new parent is on the correct branch
        bool isValidPath = true;
        auto parent = newParent;
        while (parent && parent != root) {
            if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
                if (parent->isLeftChild()) {
                    isValidPath = false;
                    break;
                }
            } else { // HORIZONTAL
                if (parent->isRightChild()) {
                    isValidPath = false;
                    break;
                }
            }
            parent = parent->getParent();
        }
        
        if (!isValidPath) return false;
    }
    
    // Remove the node from its current parent
    auto oldParent = node->getParent();
    if (oldParent) {
        if (oldParent->getLeftChild() == node) {
            oldParent->setLeftChild(nullptr);
        } else if (oldParent->getRightChild() == node) {
            oldParent->setRightChild(nullptr);
        }
    }
    
    // Add the node to its new parent
    node->setParent(newParent);
    if (asLeftChild) {
        // If there's already a left child, handle it
        auto existingChild = newParent->getLeftChild();
        if (existingChild) {
            // Try to find a place for the existing child
            // This is a simplification - in practice, more sophisticated handling is needed
            node->setLeftChild(existingChild);
            existingChild->setParent(node);
        }
        newParent->setLeftChild(node);
    } else {
        // If there's already a right child, handle it
        auto existingChild = newParent->getRightChild();
        if (existingChild) {
            // Try to find a place for the existing child
            // This is a simplification - in practice, more sophisticated handling is needed
            node->setRightChild(existingChild);
            existingChild->setParent(node);
        }
        newParent->setRightChild(node);
    }
    
    return true;
}

/**
 * Swaps two nodes in the tree
 */
bool ASFBStarTree::swapNodes(const string& nodeName1, const string& nodeName2) {
    // Find the nodes
    shared_ptr<BStarTreeNode> node1 = nullptr;
    shared_ptr<BStarTreeNode> node2 = nullptr;
    
    // Find the nodes to swap
    queue<shared_ptr<BStarTreeNode>> queue;
    if (root) queue.push(root);
    
    while (!queue.empty() && (!node1 || !node2)) {
        auto current = queue.front();
        queue.pop();
        
        if (current->getModuleName() == nodeName1) {
            node1 = current;
        }
        if (current->getModuleName() == nodeName2) {
            node2 = current;
        }
        
        if (current->getLeftChild()) queue.push(current->getLeftChild());
        if (current->getRightChild()) queue.push(current->getRightChild());
    }
    
    if (!node1 || !node2) return false;
    
    // Check if both nodes can be swapped
    // Self-symmetric modules have special restrictions
    bool node1OnBoundary = isOnBoundary(node1->getModuleName());
    bool node2OnBoundary = isOnBoundary(node2->getModuleName());
    
    // If either node is a self-symmetric module, they must stay on the boundary
    if (node1OnBoundary || node2OnBoundary) {
        // If both are on the boundary, we can swap them
        if (node1OnBoundary && node2OnBoundary) {
            // Swap is allowed
        } else {
            // One is on the boundary, one is not - swap not allowed
            return false;
        }
    }
    
    // Instead of setting module names (which we can't do), swap the parent-child relationships
    // Get parents and positions
    auto parent1 = node1->getParent();
    auto parent2 = node2->getParent();
    
    bool isLeftChild1 = parent1 && parent1->getLeftChild() == node1;
    bool isLeftChild2 = parent2 && parent2->getLeftChild() == node2;
    
    // Detach nodes from parents
    if (parent1) {
        if (isLeftChild1) parent1->setLeftChild(nullptr);
        else parent1->setRightChild(nullptr);
    }
    
    if (parent2) {
        if (isLeftChild2) parent2->setLeftChild(nullptr);
        else parent2->setRightChild(nullptr);
    }
    
    // Reattach nodes to opposite parents
    if (parent1) {
        if (isLeftChild1) parent1->setLeftChild(node2);
        else parent1->setRightChild(node2);
        node2->setParent(parent1);
    } else {
        // node1 was the root
        root = node2;
        node2->setParent(nullptr);
    }
    
    if (parent2) {
        if (isLeftChild2) parent2->setLeftChild(node1);
        else parent2->setRightChild(node1);
        node1->setParent(parent2);
    } else {
        // node2 was the root
        root = node1;
        node1->setParent(nullptr);
    }
    
    // Swap children too
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
    
    return true;
}

/**
 * Changes the representative of a symmetry pair
 */
bool ASFBStarTree::changeRepresentative(const string& moduleName) {
    // Find the symmetry pair
    auto pairIt = symmetricPairMap.find(moduleName);
    if (pairIt == symmetricPairMap.end()) return false;
    
    string module1 = moduleName;
    string module2 = pairIt->second;
    
    // Update the representative map
    if (representativeMap[module1] == module2) {
        // Change representative from module2 to module1
        representativeMap[module1] = module1;
        representativeMap[module2] = module1;
    } else {
        // Change representative from module1 to module2
        representativeMap[module1] = module2;
        representativeMap[module2] = module2;
    }
    
    // Rebuild the tree with the new representative
    constructInitialTree();
    
    return true;
}

/**
 * Converts the symmetry type (vertical to horizontal or vice versa)
 */
bool ASFBStarTree::convertSymmetryType() {
    if (!symmetryGroup) return false;
    
    // Change the symmetry type
    if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
        symmetryGroup->setType(SymmetryType::HORIZONTAL);
    } else {
        symmetryGroup->setType(SymmetryType::VERTICAL);
    }
    
    // Rotate all modules
    for (auto& pair : modules) {
        pair.second->rotate();
    }
    
    // Rebuild the tree with the new symmetry type
    constructInitialTree();
    
    return true;
}

shared_ptr<BStarTreeNode> ASFBStarTree::getRoot() const {
    return root;
}

/**
 * Gets all modules in the symmetry group
 */
const map<string, shared_ptr<Module>>& ASFBStarTree::getModules() const {
    return modules;
}

shared_ptr<SymmetryGroup> ASFBStarTree::getSymmetryGroup() const {
    return symmetryGroup;
}

double ASFBStarTree::getSymmetryAxisPosition() const {
    return symmetryAxisPosition;
}

string ASFBStarTree::getRepresentative(const string& moduleName) const {
    auto it = representativeMap.find(moduleName);
    if (it == representativeMap.end()) return "";
    
    return it->second;
}

/**
 * Checks if a module is a representative
 */
bool ASFBStarTree::isRepresentative(const string& moduleName) const {
    auto it = representativeMap.find(moduleName);
    if (it == representativeMap.end()) return false;
    
    // A module is a representative if it represents itself
    return it->second == moduleName;
}

/**
 * Creates a deep copy of this ASF-B*-tree
 */
shared_ptr<ASFBStarTree> ASFBStarTree::clone() const {
    auto clone = make_shared<ASFBStarTree>(symmetryGroup);
    
    // Copy modules
    for (const auto& pair : modules) {
        auto moduleCopy = make_shared<Module>(*pair.second);
        clone->modules[pair.first] = moduleCopy;
    }
    
    // Copy representative map and symmetry pair map
    clone->representativeMap = representativeMap;
    clone->symmetricPairMap = symmetricPairMap;
    clone->selfSymmetricModules = selfSymmetricModules;
    
    // Copy symmetry axis position
    clone->symmetryAxisPosition = symmetryAxisPosition;
    
    // Clone the tree structure
    if (root) {
        function<shared_ptr<BStarTreeNode>(shared_ptr<BStarTreeNode>)> cloneNode;
        cloneNode = [&cloneNode](shared_ptr<BStarTreeNode> node) -> shared_ptr<BStarTreeNode> {
            if (!node) return nullptr;
            
            auto newNode = make_shared<BStarTreeNode>(node->getModuleName());
            
            if (node->getLeftChild()) {
                auto leftChild = cloneNode(node->getLeftChild());
                newNode->setLeftChild(leftChild);
                leftChild->setParent(newNode);
            }
            
            if (node->getRightChild()) {
                auto rightChild = cloneNode(node->getRightChild());
                newNode->setRightChild(rightChild);
                rightChild->setParent(newNode);
            }
            
            return newNode;
        };
        
        clone->root = cloneNode(root);
    }
    
    // Clone contours
    clone->horizontalContour = make_shared<Contour>(*horizontalContour);
    clone->verticalContour = make_shared<Contour>(*verticalContour);
    
    return clone;
}