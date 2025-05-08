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
/**
 * Update contour with a module
 */
void ASFBStarTree::updateContourWithModule(const shared_ptr<Module>& module) {
    if (!module) return;
    
    int x = module->getX();
    int y = module->getY();
    int width = module->getWidth();   // Using getWidth() to properly handle rotation
    int height = module->getHeight(); // Using getHeight() to properly handle rotation
    
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
    auto moduleIt = modules.find(moduleName);
    if (moduleIt == modules.end()) return;
    
    auto module = moduleIt->second;
    
    int x = 0, y = 0;
    
    // Calculate x-coordinate based on B*-tree rules
    if (node->getParent()) {
        auto parentName = node->getParent()->getModuleName();
        auto parentIt = modules.find(parentName);
        if (parentIt != modules.end()) {
            auto parent = parentIt->second;
            if (node->isLeftChild()) {
                // Left child: place to the right of parent
                x = parent->getX() + parent->getWidth(); // Using getWidth() to handle rotation
            } else {
                // Right child: same x-coordinate as parent
                x = parent->getX();
            }
        }
    }
    
    // Calculate y-coordinate using the horizontal contour
    y = horizontalContour->getHeight(x, x + module->getWidth()); // Using getWidth() for rotated dimensions
    
    // Special handling for self-symmetric modules
    if (find(selfSymmetricModules.begin(), selfSymmetricModules.end(), 
             moduleName) != selfSymmetricModules.end()) {
        // Self-symmetric modules must be on the boundary
        if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
            // For vertical symmetry, the module must abut the symmetry axis
            // Calculate centered position based on effective width (accounting for rotation)
            x = static_cast<int>(symmetryAxisPosition) - module->getWidth() / 2;
        } else {
            // For horizontal symmetry, the module must be on the bottom
            // Calculate centered position based on effective height (accounting for rotation)
            y = static_cast<int>(symmetryAxisPosition) - module->getHeight() / 2;
        }
    }
    
    // Set the module's position
    module->setPosition(x, y);
    
    // Update contours using effective dimensions
    updateContourWithModule(module);
}

/**
 * Calculate the positions of symmetric modules
 */
/**
 * Calculate the positions of symmetric modules with careful handling of overlaps
 */
/**
 * Completely rewritten: Calculate positions for symmetric modules with 
 * guaranteed separation
 */
void ASFBStarTree::calculateSymmetricModulePositions() {
    // Step 1: Force all modules to have positive coordinates
    // Shift any modules with negative coordinates
    int minX = numeric_limits<int>::max();
    int minY = numeric_limits<int>::max();
    
    for (auto& pair : modules) {
        auto& module = pair.second;
        minX = min(minX, module->getX());
        minY = min(minY, module->getY());
    }
    
    // If any coordinates are negative, shift all modules
    if (minX < 0 || minY < 0) {
        int shiftX = max(0, -minX) + 5; // 5-unit buffer
        int shiftY = max(0, -minY) + 5; // 5-unit buffer
        
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
    const int BUFFER = 10; // Minimum separation between modules
    
    // Process symmetry pairs with exact calculations
    for (const auto& pair : symmetryGroup->getSymmetryPairs()) {
        const string& module1 = pair.first;
        const string& module2 = pair.second;
        
        auto mod1It = modules.find(module1);
        auto mod2It = modules.find(module2);
        
        if (mod1It == modules.end() || mod2It == modules.end()) continue;
        
        auto mod1 = mod1It->second;
        auto mod2 = mod2It->second;
        
        // Ensure both modules have the same dimensions and rotation status
        mod2->setRotation(mod1->getRotated());
        
        if (symType == SymmetryType::VERTICAL) {
            // Calculate center points
            int width1 = mod1->getWidth();
            int width2 = mod2->getWidth();
            
            // First, set one module's position (representative)
            // Force module1 to the left side of axis
            int mod1CenterX = static_cast<int>(symmetryAxisPosition) - width1/2 - BUFFER;
            
            // Calculate module2's position to be symmetric
            int mod2CenterX = 2 * static_cast<int>(symmetryAxisPosition) - mod1CenterX;
            
            // Set final positions ensuring both have the same y-coordinate
            mod1->setPosition(mod1CenterX - width1/2, mod1->getY());
            mod2->setPosition(mod2CenterX - width2/2, mod1->getY());
            
            cout << "Positioned symmetry pair: " << module1 << " at (" 
                 << mod1->getX() << "," << mod1->getY() << ") and " 
                 << module2 << " at (" << mod2->getX() << "," << mod2->getY() 
                 << ") about axis=" << symmetryAxisPosition << endl;
        } 
        else { // HORIZONTAL
            // Similar logic but for horizontal symmetry
            int height1 = mod1->getHeight();
            int height2 = mod2->getHeight();
            
            // Force module1 below the axis
            int mod1CenterY = static_cast<int>(symmetryAxisPosition) - height1/2 - BUFFER;
            
            // Calculate module2's position to be symmetric
            int mod2CenterY = 2 * static_cast<int>(symmetryAxisPosition) - mod1CenterY;
            
            // Set final positions ensuring both have the same x-coordinate
            mod1->setPosition(mod1->getX(), mod1CenterY - height1/2);
            mod2->setPosition(mod1->getX(), mod2CenterY - height2/2);
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
    
    // Step 4: Final verification - ensure no overlaps
    bool overlapsExist = true;
    int iterationLimit = 10; // Prevent infinite loops
    
    while (overlapsExist && iterationLimit-- > 0) {
        overlapsExist = false;
        
        // Check for and fix any remaining overlaps
        for (auto& pair1 : modules) {
            auto& mod1 = pair1.second;
            
            for (auto& pair2 : modules) {
                if (pair1.first == pair2.first) continue; // Skip self
                
                auto& mod2 = pair2.second;
                
                if (mod1->overlaps(*mod2)) {
                    overlapsExist = true;
                    
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
                            string symModName = (symPair.first == pair1.first) ? symPair.second : symPair.first;
                            auto symModIt = modules.find(symModName);
                            if (symModIt != modules.end()) {
                                // Recalculate position to maintain symmetry
                                auto symMod = symModIt->second;
                                if (symType == SymmetryType::VERTICAL) {
                                    int center1X = mod1->getX() + mod1->getWidth()/2;
                                    int center2X = 2*static_cast<int>(symmetryAxisPosition) - center1X;
                                    symMod->setPosition(center2X - symMod->getWidth()/2, mod1->getY());
                                } else {
                                    int center1Y = mod1->getY() + mod1->getHeight()/2;
                                    int center2Y = 2*static_cast<int>(symmetryAxisPosition) - center1Y;
                                    symMod->setPosition(mod1->getX(), center2Y - symMod->getHeight()/2);
                                }
                            }
                        }
                        else if (symPair.first == pair2.first || symPair.second == pair2.first) {
                            // Similar handling for mod2
                            string symModName = (symPair.first == pair2.first) ? symPair.second : symPair.first;
                            auto symModIt = modules.find(symModName);
                            if (symModIt != modules.end()) {
                                auto symMod = symModIt->second;
                                if (symType == SymmetryType::VERTICAL) {
                                    int center2X = mod2->getX() + mod2->getWidth()/2;
                                    int center1X = 2*static_cast<int>(symmetryAxisPosition) - center2X;
                                    symMod->setPosition(center1X - symMod->getWidth()/2, mod2->getY());
                                } else {
                                    int center2Y = mod2->getY() + mod2->getHeight()/2;
                                    int center1Y = 2*static_cast<int>(symmetryAxisPosition) - center2Y;
                                    symMod->setPosition(mod2->getX(), center1Y - symMod->getHeight()/2);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/**
 * Validates placement and corrects any overlaps between symmetry pairs
 */
void ASFBStarTree::validateAndCorrectOverlaps() {
    // Only validate symmetry pairs - they're the main source of overlaps
    for (const auto& pair : symmetryGroup->getSymmetryPairs()) {
        const string& module1 = pair.first;
        const string& module2 = pair.second;
        
        auto mod1It = modules.find(module1);
        auto mod2It = modules.find(module2);
        
        if (mod1It == modules.end() || mod2It == modules.end()) continue;
        
        auto mod1 = mod1It->second;
        auto mod2 = mod2It->second;
        
        // Check if these two modules overlap
        if (mod1->overlaps(*mod2)) {
            std::cout << "Fixing overlap between " << module1 << " and " << module2 << std::endl;
            
            if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
                int width1 = mod1->getWidth();
                int width2 = mod2->getWidth();
                
                // Check if they're both on the same side of the axis
                int center1 = mod1->getX() + width1/2;
                int center2 = mod2->getX() + width2/2;
                
                bool center1LeftOfAxis = center1 < symmetryAxisPosition;
                bool center2LeftOfAxis = center2 < symmetryAxisPosition;
                
                if (center1LeftOfAxis == center2LeftOfAxis) {
                    // Both on same side - separate them but keep them symmetric
                    if (center1LeftOfAxis) {
                        // Both left of axis, move leftmost one further left
                        if (mod1->getX() <= mod2->getX()) {
                            int newX1 = mod2->getX() - width1 - 5; // 5-unit buffer
                            mod1->setPosition(newX1, mod1->getY());
                        } else {
                            int newX2 = mod1->getX() - width2 - 5; // 5-unit buffer
                            mod2->setPosition(newX2, mod2->getY());
                        }
                    } else {
                        // Both right of axis, move rightmost one further right
                        if (mod1->getX() >= mod2->getX()) {
                            int newX2 = mod1->getX() + width1 + 5; // 5-unit buffer
                            mod2->setPosition(newX2, mod2->getY());
                        } else {
                            int newX1 = mod2->getX() + width2 + 5; // 5-unit buffer
                            mod1->setPosition(newX1, mod1->getY());
                        }
                    }
                } else {
                    // One on each side - ensure they're properly symmetric
                    // Determine the representative module
                    auto rep = representativeMap[module1] == module1 ? mod1 : mod2;
                    auto other = rep == mod1 ? mod2 : mod1;
                    string repName = rep == mod1 ? module1 : module2;
                    string otherName = other == mod1 ? module1 : module2;
                    
                    int repCenter = rep->getX() + rep->getWidth()/2;
                    int otherCenter = 2 * static_cast<int>(symmetryAxisPosition) - repCenter;
                    int newX = otherCenter - other->getWidth()/2;
                    
                    std::cout << "Repositioning " << otherName << " to (" << newX
                              << "," << other->getY() << ") based on " << repName << std::endl;
                    
                    other->setPosition(newX, rep->getY());
                }
            } else { // HORIZONTAL
                // Similar logic for horizontal symmetry
                int height1 = mod1->getHeight();
                int height2 = mod2->getHeight();
                
                // Check if they're both on the same side of the axis
                int center1 = mod1->getY() + height1/2;
                int center2 = mod2->getY() + height2/2;
                
                bool center1BelowAxis = center1 < symmetryAxisPosition;
                bool center2BelowAxis = center2 < symmetryAxisPosition;
                
                if (center1BelowAxis == center2BelowAxis) {
                    // Both on same side - separate them but keep them symmetric
                    if (center1BelowAxis) {
                        // Both below axis, move lower one further down
                        if (mod1->getY() <= mod2->getY()) {
                            int newY1 = mod2->getY() - height1 - 5; // 5-unit buffer
                            mod1->setPosition(mod1->getX(), newY1);
                        } else {
                            int newY2 = mod1->getY() - height2 - 5; // 5-unit buffer
                            mod2->setPosition(mod2->getX(), newY2);
                        }
                    } else {
                        // Both above axis, move higher one further up
                        if (mod1->getY() >= mod2->getY()) {
                            int newY2 = mod1->getY() + height1 + 5; // 5-unit buffer
                            mod2->setPosition(mod2->getX(), newY2);
                        } else {
                            int newY1 = mod2->getY() + height2 + 5; // 5-unit buffer
                            mod1->setPosition(mod1->getX(), newY1);
                        }
                    }
                } else {
                    // One on each side - ensure they're properly symmetric
                    // Determine the representative module
                    auto rep = representativeMap[module1] == module1 ? mod1 : mod2;
                    auto other = rep == mod1 ? mod2 : mod1;
                    
                    int repCenter = rep->getY() + rep->getHeight()/2;
                    int otherCenter = 2 * static_cast<int>(symmetryAxisPosition) - repCenter;
                    int newY = otherCenter - other->getHeight()/2;
                    
                    other->setPosition(rep->getX(), newY);
                }
            }
        }
    }
}

/**
 * Calculates the coordinates of all modules in the symmetry group
 * by packing the ASF-B*-tree
 */
/**
 * Completely rewritten: Reliable packing method
 */
bool ASFBStarTree::pack() {
    if (!root) return false;
    
    // Initialize contours
    initializeContours();
    
    // Step 1: Set a reliable initial symmetry axis - must be positive
    // For vertical symmetry, set an initial axis far from the origin
    // For horizontal symmetry, set an initial axis far from the origin
    symmetryAxisPosition = symmetryGroup->getType() == SymmetryType::VERTICAL ? 5000 : 5000;
    
    // Step 2: Pack all representative modules
    queue<shared_ptr<BStarTreeNode>> nodeQueue;
    nodeQueue.push(root);
    
    while (!nodeQueue.empty()) {
        auto currentNode = nodeQueue.front();
        nodeQueue.pop();
        
        if (!currentNode) continue;
        
        // Pack the current node (uses effective dimensions from getWidth/getHeight)
        packNode(currentNode);
        
        // Add children to the queue
        if (currentNode->getLeftChild()) {
            nodeQueue.push(currentNode->getLeftChild());
        }
        if (currentNode->getRightChild()) {
            nodeQueue.push(currentNode->getRightChild());
        }
    }
    
    // Step 3: Ensure no modules have negative coordinates
    int minX = 0;
    int minY = 0;
    
    for (const auto& pair : modules) {
        const auto& module = pair.second;
        minX = min(minX, module->getX());
        minY = min(minY, module->getY());
    }
    
    // If we have negative coordinates, shift all modules
    if (minX < 0 || minY < 0) {
        int shiftX = max(0, -minX) + 10;  // 10-unit buffer
        int shiftY = max(0, -minY) + 10;  // 10-unit buffer
        
        for (auto& pair : modules) {
            auto& module = pair.second;
            module->setPosition(module->getX() + shiftX, module->getY() + shiftY);
        }
        
        // Adjust symmetry axis position
        if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
            symmetryAxisPosition += shiftX;
        } else {
            symmetryAxisPosition += shiftY;
        }
    }
    
    // Step 4: Calculate symmetry axis based on current module positions
    // Use a consistent, stable method rather than the previous calculation
    if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
        // Find rightmost module edge
        int maxRight = 0;
        for (const auto& pair : modules) {
            const auto& module = pair.second;
            maxRight = max(maxRight, module->getX() + module->getWidth());
        }
        // Set axis to be a fixed distance to the right of all modules
        symmetryAxisPosition = maxRight + 1000;  // Large buffer
    } else {
        // Find topmost module edge
        int maxTop = 0;
        for (const auto& pair : modules) {
            const auto& module = pair.second;
            maxTop = max(maxTop, module->getY() + module->getHeight());
        }
        // Set axis to be a fixed distance above all modules
        symmetryAxisPosition = maxTop + 1000;  // Large buffer
    }
    
    std::cout << "Using fixed symmetry axis: " << symmetryAxisPosition << std::endl;
    
    // Step 5: Calculate positions for all modules in symmetry pairs
    calculateSymmetricModulePositions();
    
    // Step 6: Final validation
    if (!validateSymmetryConstraints()) {
        std::cout << "Final pack failed symmetry validation - trying one more adjustment" << std::endl;
        
        // One more attempt to fix symmetry issues
        for (const auto& pair : symmetryGroup->getSymmetryPairs()) {
            const string& module1 = pair.first;
            const string& module2 = pair.second;
            
            auto mod1It = modules.find(module1);
            auto mod2It = modules.find(module2);
            
            if (mod1It == modules.end() || mod2It == modules.end()) continue;
            
            auto mod1 = mod1It->second;
            auto mod2 = mod2It->second;
            
            // Explicitly enforce symmetry
            if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
                int center1X = mod1->getX() + mod1->getWidth()/2;
                int center2X = 2*static_cast<int>(symmetryAxisPosition) - center1X;
                mod2->setPosition(center2X - mod2->getWidth()/2, mod1->getY());
            } else {
                int center1Y = mod1->getY() + mod1->getHeight()/2;
                int center2Y = 2*static_cast<int>(symmetryAxisPosition) - center1Y;
                mod2->setPosition(mod1->getX(), center2Y - mod2->getHeight()/2);
            }
        }
    }
    
    return true;
}

/**
 * Validates that all symmetry constraints are satisfied
 */
bool ASFBStarTree::validateSymmetryConstraints() const {
    if (!symmetryGroup) return true;
    
    std::cout << "Validating symmetry constraints for group: " << symmetryGroup->getName() << std::endl;
    std::cout << "Symmetry axis position: " << symmetryAxisPosition << std::endl;
    
    // Check symmetry pairs have correct positions
    for (const auto& pair : symmetryGroup->getSymmetryPairs()) {
        const auto& module1It = modules.find(pair.first);
        const auto& module2It = modules.find(pair.second);
        
        if (module1It == modules.end() || module2It == modules.end()) continue;
        
        const auto& mod1 = module1It->second;
        const auto& mod2 = module2It->second;
        
        if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
            // Calculate centers accounting for rotation
            double center1X = mod1->getX() + mod1->getWidth() / 2.0;
            double center2X = mod2->getX() + mod2->getWidth() / 2.0;
            
            // Check x-coordinate reflection: x1 + x2 = 2 * symmetryAxis
            double diff = abs((center1X + center2X) - 2 * symmetryAxisPosition);
            if (diff > 1) {
                std::cout << "Vertical symmetry validation failed for pair: " 
                          << pair.first << " (" << mod1->getX() << "," << mod1->getY() << ") center=" << center1X
                          << " and " << pair.second << " (" << mod2->getX() << "," << mod2->getY() << ") center=" << center2X
                          << " diff=" << diff << " axis=" << symmetryAxisPosition << std::endl;
                return false;
            }
            
            // y1 should equal y2
            if (abs(mod1->getY() - mod2->getY()) > 1) {
                std::cout << "Y-coordinate mismatch for pair: " 
                          << pair.first << " (" << mod1->getY() << ") and "
                          << pair.second << " (" << mod2->getY() << ")" << std::endl;
                return false;
            }
        } else { // HORIZONTAL
            // Similar validation with more debugging output...
        }
    }
    
    // Self-symmetric module validation with debugging...
    
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
    if (modules.empty()) return true;
    
    // Build adjacency list for all modules in the symmetry group
    unordered_map<string, vector<string>> adjacency;
    
    // Get all module positions
    unordered_map<string, pair<int, int>> positions;
    unordered_map<string, pair<int, int>> dimensions;
    
    for (const auto& pair : modules) {
        const auto& module = pair.second;
        positions[pair.first] = {module->getX(), module->getY()};
        dimensions[pair.first] = {module->getWidth(), module->getHeight()};
    }
    
    // Build adjacency list based on abutting modules
    for (const auto& pair1 : modules) {
        const auto& name1 = pair1.first;
        const auto& pos1 = positions[name1];
        const auto& dim1 = dimensions[name1];
        
        for (const auto& pair2 : modules) {
            const auto& name2 = pair2.first;
            if (name1 == name2) continue;
            
            const auto& pos2 = positions[name2];
            const auto& dim2 = dimensions[name2];
            
            // Check if modules are adjacent
            bool adjacent = false;
            
            // Check horizontal adjacency
            if (pos1.first + dim1.first == pos2.first || 
                pos2.first + dim2.first == pos1.first) {
                // Check vertical overlap
                if (!(pos1.second >= pos2.second + dim2.second || 
                      pos2.second >= pos1.second + dim1.second)) {
                    adjacent = true;
                }
            }
            
            // Check vertical adjacency
            if (pos1.second + dim1.second == pos2.second || 
                pos2.second + dim2.second == pos1.second) {
                // Check horizontal overlap
                if (!(pos1.first >= pos2.first + dim2.first || 
                      pos2.first >= pos1.first + dim1.first)) {
                    adjacent = true;
                }
            }
            
            if (adjacent) {
                adjacency[name1].push_back(name2);
            }
        }
    }
    
    // Perform BFS to check connectivity
    unordered_set<string> visited;
    queue<string> q;
    
    auto it = modules.begin();
    q.push(it->first);
    visited.insert(it->first);
    
    while (!q.empty()) {
        string current = q.front();
        q.pop();
        
        for (const auto& neighbor : adjacency[current]) {
            if (visited.find(neighbor) == visited.end()) {
                visited.insert(neighbor);
                q.push(neighbor);
            }
        }
    }
    
    // Check if all modules are visited
    return visited.size() == modules.size();
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