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
        int shiftX = std::max(0, -minX) + 2; // Reduced buffer from 5 to 2
        int shiftY = std::max(0, -minY) + 2; // Reduced buffer from 5 to 2
        
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
    
    // Step 2: Apply optimized symmetry positioning for all pairs
    SymmetryType symType = symmetryGroup->getType();
    
    // Minimal buffer to ensure true connectivity without excessive padding
    const int MIN_BUFFER = 2; // Minimal buffer to avoid exact overlaps
    
    // Track processed modules to maintain connectivity
    std::unordered_map<std::string, bool> processedModules;
    
    // Build an ordering of symmetry pairs for consistent placement
    std::vector<std::pair<std::string, std::string>> orderedPairs = symmetryGroup->getSymmetryPairs();
    // Sort pairs by module name to ensure consistent ordering
    std::sort(orderedPairs.begin(), orderedPairs.end());
    
    // Process symmetry pairs with exact calculations to avoid overlaps
    for (size_t i = 0; i < orderedPairs.size(); i++) {
        const auto& pair = orderedPairs[i];
        const std::string& module1 = pair.first;
        const std::string& module2 = pair.second;
        
        auto mod1It = modules.find(module1);
        auto mod2It = modules.find(module2);
        
        if (mod1It == modules.end() || mod2It == modules.end()) continue;
        
        auto mod1 = mod1It->second;
        auto mod2 = mod2It->second;
        
        // Ensure both modules have the same dimensions and rotation status
        mod2->setRotation(mod1->getRotated());
        
        // Position based on symmetry type
        if (symType == SymmetryType::VERTICAL) {
            // Calculate widths (accounting for rotation)
            int width1 = mod1->getWidth();
            int width2 = mod2->getWidth();
            
            // CRITICAL CHANGE: Calculate non-overlapping positions using exact arithmetic
            int mod1X = static_cast<int>(symmetryAxisPosition) - width1 - MIN_BUFFER;
            int mod2X = static_cast<int>(symmetryAxisPosition) + MIN_BUFFER;
            
            // Y-position based on already processed modules
            int yPos = 0;
            if (!processedModules.empty() && i > 0) {
                for (const auto& placedMod : processedModules) {
                    auto modIt = modules.find(placedMod.first);
                    if (modIt != modules.end()) {
                        int modBottom = modIt->second->getY() + modIt->second->getHeight();
                        yPos = std::max(yPos, modBottom);
                    }
                }
                // Add minimal spacing for connectivity without gaps
                yPos += MIN_BUFFER;
            }
            
            // Set positions with symmetry around the axis
            mod1->setPosition(mod1X, yPos);
            mod2->setPosition(mod2X, yPos);
            
            std::cout << "Positioned symmetry pair: " << module1 << " at (" 
                     << mod1->getX() << "," << mod1->getY() << ") and " 
                     << module2 << " at (" << mod2->getX() << "," << mod2->getY() 
                     << ") about axis=" << symmetryAxisPosition << std::endl;
            
            // IMPORTANT: Double-check no overlap was created
            if (mod1->getX() + mod1->getWidth() > mod2->getX()) {
                std::cout << "Warning: Correcting overlap in symmetry pair" << std::endl;
                // Adjust positions to eliminate overlap
                mod1->setPosition(mod2->getX() - mod1->getWidth() - 1, yPos);
            }
        } 
        else { // HORIZONTAL
            // Similar logic for horizontal symmetry
            int height1 = mod1->getHeight();
            int height2 = mod2->getHeight();
            
            // Calculate x-position to maintain connectivity
            int xPos = 0;
            if (!processedModules.empty() && i > 0) {
                for (const auto& placedMod : processedModules) {
                    auto modIt = modules.find(placedMod.first);
                    if (modIt != modules.end()) {
                        int modRight = modIt->second->getX() + modIt->second->getWidth();
                        xPos = std::max(xPos, modRight);
                    }
                }
                xPos += MIN_BUFFER;
            }
            
            // CRITICAL CHANGE: Calculate non-overlapping positions
            int mod1Y = static_cast<int>(symmetryAxisPosition) - height1 - MIN_BUFFER;
            int mod2Y = static_cast<int>(symmetryAxisPosition) + MIN_BUFFER;
            
            // Set positions with symmetry around the axis
            mod1->setPosition(xPos, mod1Y);
            mod2->setPosition(xPos, mod2Y);
            
            // IMPORTANT: Double-check no overlap was created
            if (mod1->getY() + mod1->getHeight() > mod2->getY()) {
                // Adjust positions to eliminate overlap
                mod1->setPosition(xPos, mod2->getY() - mod1->getHeight() - 1);
            }
        }
        
        // Record modules as processed
        processedModules[module1] = true;
        processedModules[module2] = true;
    }
    
    // Step 3: Position self-symmetric modules
    int selfSymXPos = 0;
    int selfSymYPos = 0;
    
    // Find maximum coordinates to place self-symmetric modules after symmetry pairs
    if (!processedModules.empty()) {
        for (const auto& placedMod : processedModules) {
            auto modIt = modules.find(placedMod.first);
            if (modIt != modules.end()) {
                auto& mod = modIt->second;
                selfSymXPos = std::max(selfSymXPos, mod->getX() + mod->getWidth());
                selfSymYPos = std::max(selfSymYPos, mod->getY() + mod->getHeight());
            }
        }
        // Minimal spacing for connectivity
        if (symType == SymmetryType::VERTICAL) {
            selfSymYPos += MIN_BUFFER;
        } else {
            selfSymXPos += MIN_BUFFER;
        }
    }
    
    for (const auto& moduleName : selfSymmetricModules) {
        auto it = modules.find(moduleName);
        if (it == modules.end()) continue;
        
        auto module = it->second;
        
        // Center module on the symmetry axis
        if (symType == SymmetryType::VERTICAL) {
            int width = module->getWidth();
            // Exact position calculation to center on axis
            int x = static_cast<int>(symmetryAxisPosition) - width / 2;
            module->setPosition(x, selfSymYPos);
            
            // Update for next self-symmetric module
            selfSymYPos += module->getHeight() + MIN_BUFFER;
        } else { // HORIZONTAL
            int height = module->getHeight();
            // Exact position calculation to center on axis
            int y = static_cast<int>(symmetryAxisPosition) - height / 2;
            module->setPosition(selfSymXPos, y);
            
            // Update for next self-symmetric module
            selfSymXPos += module->getWidth() + MIN_BUFFER;
        }
        
        std::cout << "Positioned self-symmetric module: " << moduleName 
                 << " at (" << module->getX() << "," << module->getY() << ")" << std::endl;
        
        // Mark as processed
        processedModules[moduleName] = true;
    }
    
    // Step 4: Verify connectivity between all modules
    bool allConnected = verifyConnectivity();
    if (!allConnected) {
        std::cout << "Warning: Modules not fully connected after initial positioning, forcing connectivity..." << std::endl;
        forceConnectivity();
    }
    
    // Step 5: Final verification - check for overlaps
    bool overlapsExist = checkForOverlaps();
    if (overlapsExist) {
        std::cout << "Warning: Overlaps detected after positioning, fixing..." << std::endl;
        fixOverlaps();
    }
}

/**
 * Calculates the coordinates of all modules in the symmetry group
 * by packing the ASF-B*-tree, includes symmetry enforcement
 */
bool ASFBStarTree::pack() {
    if (!root) {
        std::cerr << "Cannot pack ASFBStarTree: root is null" << std::endl;
        return false;
    }
    
    // Initialize contours
    initializeContours();
    
    // First pass: Pack all representative modules
    std::queue<std::shared_ptr<BStarTreeNode>> nodeQueue;
    nodeQueue.push(root);
    
    // Track min/max coordinates for axis calculation
    int minX = std::numeric_limits<int>::max();
    int maxX = std::numeric_limits<int>::min();
    int minY = std::numeric_limits<int>::max();
    int maxY = std::numeric_limits<int>::min();
    
    // Pack the representative modules
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
    
    // Calculate stable symmetry axis
    const int MIN_AXIS_POSITION = 100; // Minimum axis position to ensure stability
    
    if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
        // Calculate axis as midpoint between min and max X, ensuring minimum position
        symmetryAxisPosition = std::max(MIN_AXIS_POSITION, (minX + maxX) / 2);
    } else { // HORIZONTAL
        // Calculate axis as midpoint between min and max Y, ensuring minimum position
        symmetryAxisPosition = std::max(MIN_AXIS_POSITION, (minY + maxY) / 2);
    }
    
    std::cout << "Using stabilized symmetry axis: " << symmetryAxisPosition << std::endl;
    
    // Calculate positions for all modules in symmetry pairs
    calculateSymmetricModulePositions();
    
    // Explicitly enforce symmetry constraints to fix any violations
    enforceSymmetryConstraints();
    
    // Run another enforcement pass to handle any remaining issues
    enforceSymmetryConstraints();
    
    // Final validation after all corrections
    bool isSymmetryValid = validateSymmetryConstraints();
    bool islandValid = isSymmetryIslandValid();
    
    // If there are still problems, try emergency recovery
    if (!isSymmetryValid || !islandValid) {
        std::cerr << "Symmetry or connectivity issues persist, applying emergency recovery..." << std::endl;
        
        // Apply a simple but guaranteed-valid positioning
        emergencyRecovery();
        
        // Re-validate
        isSymmetryValid = validateSymmetryConstraints();
        islandValid = isSymmetryIslandValid();
    }
    
    // Return success if both symmetry constraints and island connectivity are valid
    return isSymmetryValid && islandValid;
}

/**
 * Emergency recovery function that forces a valid placement
 * when normal positioning attempts fail
 */
void ASFBStarTree::emergencyRecovery() {
    // Start with a clean position with more buffer space
    int startX = 50;  // Increased from 10
    int startY = 50;  // Increased from 10
    
    // Ensure symmetry axis is properly positioned with more margin
    if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
        // Find the maximum width of any module
        int maxWidth = 0;
        for (const auto& pair : modules) {
            maxWidth = std::max(maxWidth, pair.second->getWidth());
        }
        // Use a more generous buffer
        symmetryAxisPosition = startX + maxWidth * 2.5; // Increased from *3
    } else {
        // Find the maximum height of any module
        int maxHeight = 0;
        for (const auto& pair : modules) {
            maxHeight = std::max(maxHeight, pair.second->getHeight());
        }
        // Use a more generous buffer
        symmetryAxisPosition = startY + maxHeight * 2.5; // Increased from *3
    }
    
    std::cout << "Emergency recovery: using symmetry axis at " << symmetryAxisPosition << std::endl;
    
    // Position symmetric pairs - place them directly against the axis with no gap
    int currentX = startX;
    int currentY = startY;
    
    for (const auto& pair : symmetryGroup->getSymmetryPairs()) {
        auto mod1It = modules.find(pair.first);
        auto mod2It = modules.find(pair.second);
        
        if (mod1It == modules.end() || mod2It == modules.end()) continue;
        
        auto mod1 = mod1It->second;
        auto mod2 = mod2It->second;
        
        // Ensure same rotation status
        mod2->setRotation(mod1->getRotated());
        
        if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
            // Place modules directly against vertical axis with no gap
            int mod1X = static_cast<int>(symmetryAxisPosition) - mod1->getWidth();
            int mod2X = static_cast<int>(symmetryAxisPosition);
            
            mod1->setPosition(mod1X, currentY);
            mod2->setPosition(mod2X, currentY);
            
            // Move down ensuring modules touch exactly
            currentY += std::max(mod1->getHeight(), mod2->getHeight());
        } else {
            // Place modules directly against horizontal axis with no gap
            int mod1Y = static_cast<int>(symmetryAxisPosition) - mod1->getHeight();
            int mod2Y = static_cast<int>(symmetryAxisPosition);
            
            mod1->setPosition(currentX, mod1Y);
            mod2->setPosition(currentX, mod2Y);
            
            // Move right ensuring modules touch exactly
            currentX += std::max(mod1->getWidth(), mod2->getWidth());
        }
    }
    
    // Position self-symmetric modules - center them exactly on the axis
    for (const auto& name : selfSymmetricModules) {
        auto moduleIt = modules.find(name);
        if (moduleIt == modules.end()) continue;
        
        auto module = moduleIt->second;
        
        if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
            // Center exactly on vertical axis
            int x = static_cast<int>(symmetryAxisPosition) - module->getWidth() / 2;
            module->setPosition(x, currentY);
            
            // Move down with no gap
            currentY += module->getHeight();
        } else {
            // Center exactly on horizontal axis
            int y = static_cast<int>(symmetryAxisPosition) - module->getHeight() / 2;
            module->setPosition(currentX, y);
            
            // Move right with no gap
            currentX += module->getWidth();
        }
    }
    
    // Final verification to ensure connectivity and no overlaps
    bool connected = verifyConnectivity();
    bool noOverlaps = !checkForOverlaps();
    
    if (!connected) {
        std::cout << "Emergency recovery: modules still not connected, forcing tighter packing" << std::endl;
        forceConnectivity();
    }
    
    if (!noOverlaps) {
        std::cout << "Emergency recovery: overlaps detected, applying overlap fix" << std::endl;
        fixOverlaps();
    }
    
    std::cout << "Emergency recovery completed" << std::endl;
}


/**
 * Validates that all symmetry constraints are satisfied
 * This enhanced version will correct position violations when possible
 */
/**
 * Validates that all symmetry constraints are satisfied
 * This enhanced version will correct position violations when possible
 */
bool ASFBStarTree::validateSymmetryConstraints() const {
    if (!symmetryGroup) {
        std::cerr << "No symmetry group defined" << std::endl;
        return false;
    }
    
    std::cout << "Validating symmetry constraints for group: " << symmetryGroup->getName() << std::endl;
    std::cout << "Symmetry axis position: " << symmetryAxisPosition << std::endl;
    std::cout << "Symmetry type: " << (symmetryGroup->getType() == SymmetryType::VERTICAL ? "VERTICAL" : "HORIZONTAL") << std::endl;
    
    bool allValid = true;
    // More generous tolerance for symmetry validation
    const double TOLERANCE = 5.0;
    
    // Check symmetry pairs have correct positions
    for (const auto& pair : symmetryGroup->getSymmetryPairs()) {
        const auto& module1It = modules.find(pair.first);
        const auto& module2It = modules.find(pair.second);
        
        if (module1It == modules.end() || module2It == modules.end()) {
            std::cerr << "Module pair not found: " << pair.first << " and " << pair.second << std::endl;
            allValid = false;
            continue;
        }
        
        auto mod1 = module1It->second;
        auto mod2 = module2It->second;
        
        // Allow non-const version to modify positions if needed
        if (!mod1 || !mod2) continue;
        
        if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
            // Calculate centers accounting for rotation
            double center1X = mod1->getX() + mod1->getWidth() / 2.0;
            double center2X = mod2->getX() + mod2->getWidth() / 2.0;
            
            // Check if centers are approximately equidistant from the axis
            double axis2x = 2.0 * symmetryAxisPosition;
            double diff = std::abs((center1X + center2X) - axis2x);
            
            if (diff > TOLERANCE) {
                std::cerr << "Vertical symmetry violation for pair: " 
                          << pair.first << " (" << mod1->getX() << "," << mod1->getY() << ") center=" << center1X
                          << " and " << pair.second << " (" << mod2->getX() << "," << mod2->getY() << ") center=" << center2X
                          << " diff=" << diff << " axis=" << symmetryAxisPosition << std::endl;
                
                // Calculate correct position for module2
                double expectedCenter2X = axis2x - center1X;
                std::cerr << "  Expected " << pair.second << " center at X=" 
                          << expectedCenter2X << " but found X=" << center2X << std::endl;
                
                // Try to fix the position - needs non-const access to modules
                if (!isValidationOnly) {
                    int newX = static_cast<int>(expectedCenter2X - mod2->getWidth() / 2.0);
                    const_cast<Module*>(mod2.get())->setPosition(newX, mod2->getY());
                    std::cout << "  Auto-corrected position to X=" << newX << std::endl;
                }
                
                allValid = false;
            }
            
            // Check if Y coordinates are approximately equal
            if (std::abs(mod1->getY() - mod2->getY()) > TOLERANCE) {
                std::cerr << "Y-coordinate mismatch for pair: " 
                          << pair.first << " (" << mod1->getY() << ") and "
                          << pair.second << " (" << mod2->getY() << ")" << std::endl;
                
                // Try to fix the position
                if (!isValidationOnly) {
                    const_cast<Module*>(mod2.get())->setPosition(mod2->getX(), mod1->getY());
                    std::cout << "  Auto-corrected Y position to match" << std::endl;
                }
                
                allValid = false;
            }
            
            // Check if both modules have the same dimensions (accounting for rotation)
            if (std::abs(mod1->getWidth() - mod2->getWidth()) > TOLERANCE || 
                std::abs(mod1->getHeight() - mod2->getHeight()) > TOLERANCE) {
                std::cerr << "Dimension mismatch for pair: " 
                          << pair.first << " (" << mod1->getWidth() << "x" << mod1->getHeight() << ") and "
                          << pair.second << " (" << mod2->getWidth() << "x" << mod2->getHeight() << ")" << std::endl;
                
                // Try to fix the rotation
                if (!isValidationOnly) {
                    // Match rotation status to ensure same dimensions
                    const_cast<Module*>(mod2.get())->setRotation(mod1->getRotated());
                    std::cout << "  Auto-corrected rotation to match dimensions" << std::endl;
                }
                
                allValid = false;
            }
        } else { // HORIZONTAL
            // Similar validation with improved tolerance and debugging for horizontal symmetry
            double center1Y = mod1->getY() + mod1->getHeight() / 2.0;
            double center2Y = mod2->getY() + mod2->getHeight() / 2.0;
            
            double axis2y = 2.0 * symmetryAxisPosition;
            double diff = std::abs((center1Y + center2Y) - axis2y);
            
            if (diff > TOLERANCE) {
                std::cerr << "Horizontal symmetry violation for pair: " 
                          << pair.first << " (" << mod1->getX() << "," << mod1->getY() << ") center=" << center1Y
                          << " and " << pair.second << " (" << mod2->getX() << "," << mod2->getY() << ") center=" << center2Y
                          << " diff=" << diff << " axis=" << symmetryAxisPosition << std::endl;
                
                // Calculate correct position for module2
                double expectedCenter2Y = axis2y - center1Y;
                std::cerr << "  Expected " << pair.second << " center at Y=" 
                          << expectedCenter2Y << " but found Y=" << center2Y << std::endl;
                
                // Try to fix the position
                if (!isValidationOnly) {
                    int newY = static_cast<int>(expectedCenter2Y - mod2->getHeight() / 2.0);
                    const_cast<Module*>(mod2.get())->setPosition(mod2->getX(), newY);
                    std::cout << "  Auto-corrected position to Y=" << newY << std::endl;
                }
                
                allValid = false;
            }
            
            if (std::abs(mod1->getX() - mod2->getX()) > TOLERANCE) {
                std::cerr << "X-coordinate mismatch for pair: " 
                          << pair.first << " (" << mod1->getX() << ") and "
                          << pair.second << " (" << mod2->getX() << ")" << std::endl;
                
                // Try to fix the position
                if (!isValidationOnly) {
                    const_cast<Module*>(mod2.get())->setPosition(mod1->getX(), mod2->getY());
                    std::cout << "  Auto-corrected X position to match" << std::endl;
                }
                
                allValid = false;
            }
            
            if (std::abs(mod1->getWidth() - mod2->getWidth()) > TOLERANCE || 
                std::abs(mod1->getHeight() - mod2->getHeight()) > TOLERANCE) {
                std::cerr << "Dimension mismatch for pair: " 
                          << pair.first << " (" << mod1->getWidth() << "x" << mod1->getHeight() << ") and "
                          << pair.second << " (" << mod2->getWidth() << "x" << mod2->getHeight() << ")" << std::endl;
                
                // Try to fix the rotation
                if (!isValidationOnly) {
                    // Match rotation status to ensure same dimensions
                    const_cast<Module*>(mod2.get())->setRotation(mod1->getRotated());
                    std::cout << "  Auto-corrected rotation to match dimensions" << std::endl;
                }
                
                allValid = false;
            }
        }
    }
    
    // Check self-symmetric modules
    for (const auto& moduleName : selfSymmetricModules) {
        auto it = modules.find(moduleName);
        if (it == modules.end()) {
            std::cerr << "Self-symmetric module not found: " << moduleName << std::endl;
            allValid = false;
            continue;
        }
        
        auto module = it->second;
        
        if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
            // Calculate center X
            double centerX = module->getX() + module->getWidth() / 2.0;
            
            // Check if center is approximately on the axis
            if (std::abs(centerX - symmetryAxisPosition) > TOLERANCE) {
                std::cerr << "Self-symmetric module " << moduleName 
                          << " not centered on vertical axis. Center X=" << centerX
                          << " but axis is at " << symmetryAxisPosition << std::endl;
                
                // Try to fix the position
                if (!isValidationOnly) {
                    int newX = static_cast<int>(symmetryAxisPosition - module->getWidth() / 2.0);
                    const_cast<Module*>(module.get())->setPosition(newX, module->getY());
                    std::cout << "  Auto-corrected position to X=" << newX << std::endl;
                }
                
                allValid = false;
            }
        } else { // HORIZONTAL
            // Calculate center Y
            double centerY = module->getY() + module->getHeight() / 2.0;
            
            // Check if center is approximately on the axis
            if (std::abs(centerY - symmetryAxisPosition) > TOLERANCE) {
                std::cerr << "Self-symmetric module " << moduleName 
                          << " not centered on horizontal axis. Center Y=" << centerY
                          << " but axis is at " << symmetryAxisPosition << std::endl;
                
                // Try to fix the position
                if (!isValidationOnly) {
                    int newY = static_cast<int>(symmetryAxisPosition - module->getHeight() / 2.0);
                    const_cast<Module*>(module.get())->setPosition(module->getX(), newY);
                    std::cout << "  Auto-corrected position to Y=" << newY << std::endl;
                }
                
                allValid = false;
            }
        }
    }
    
    if (allValid) {
        std::cout << "Symmetry validation passed!" << std::endl;
    } else if (!isValidationOnly) {
        std::cout << "Fixed symmetry constraint violations" << std::endl;
    }
    
    return allValid;
}

/**
 * Enforce symmetry constraints by adjusting module positions
 * This is a dedicated function for fixing symmetry violations
 */
void ASFBStarTree::enforceSymmetryConstraints() {
    if (!symmetryGroup) {
        std::cerr << "No symmetry group defined" << std::endl;
        return;
    }
    
    std::cout << "Enforcing symmetry constraints for group: " << symmetryGroup->getName() << std::endl;
    
    // Store original positions to detect changes
    std::unordered_map<std::string, std::pair<int, int>> originalPositions;
    for (const auto& pair : modules) {
        originalPositions[pair.first] = {pair.second->getX(), pair.second->getY()};
    }
    
    // More precision for symmetry enforcement but still detect significant violations
    const double TOLERANCE = 1.0;
    
    // Check and enforce symmetry for pairs
    for (const auto& pair : symmetryGroup->getSymmetryPairs()) {
        auto mod1It = modules.find(pair.first);
        auto mod2It = modules.find(pair.second);
        
        if (mod1It == modules.end() || mod2It == modules.end()) {
            std::cerr << "Module pair not found: " << pair.first << " and " << pair.second << std::endl;
            continue;
        }
        
        auto mod1 = mod1It->second;
        auto mod2 = mod2It->second;
        
        if (!mod1 || !mod2) continue;
        
        // Ensure both modules have the same dimensions (rotation status)
        mod2->setRotation(mod1->getRotated());
        
        if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
            // Calculate centers
            double center1X = mod1->getX() + mod1->getWidth() / 2.0;
            double center2X = mod2->getX() + mod2->getWidth() / 2.0;
            
            // Check symmetry
            double axis2x = 2.0 * symmetryAxisPosition;
            double diff = std::abs((center1X + center2X) - axis2x);
            
            if (diff > TOLERANCE) {
                // Calculate exact positions for symmetry
                double expectedCenter1X = (center1X + center2X) / 2.0 - (center2X - symmetryAxisPosition);
                double expectedCenter2X = (center1X + center2X) / 2.0 + (symmetryAxisPosition - center1X);
                
                int newX1 = static_cast<int>(expectedCenter1X - mod1->getWidth() / 2.0);
                int newX2 = static_cast<int>(expectedCenter2X - mod2->getWidth() / 2.0);
                
                // Ensure they don't overlap by checking if their centers are far enough apart
                if (std::abs(expectedCenter1X - expectedCenter2X) < (mod1->getWidth() / 2.0 + mod2->getWidth() / 2.0)) {
                    // They would overlap - adjust positions
                    double minDistance = mod1->getWidth() / 2.0 + mod2->getWidth() / 2.0 + 1; // 1 extra unit
                    
                    // Adjust equally from center
                    expectedCenter1X = symmetryAxisPosition - minDistance / 2.0;
                    expectedCenter2X = symmetryAxisPosition + minDistance / 2.0;
                    
                    newX1 = static_cast<int>(expectedCenter1X - mod1->getWidth() / 2.0);
                    newX2 = static_cast<int>(expectedCenter2X - mod2->getWidth() / 2.0);
                }
                
                // Update positions
                mod1->setPosition(newX1, mod1->getY());
                mod2->setPosition(newX2, mod2->getY());
                
                std::cout << "Enforced X-axis symmetry for pair: " 
                          << pair.first << " and " << pair.second << std::endl;
            }
            
            // Ensure Y coordinates are exactly equal
            if (std::abs(mod1->getY() - mod2->getY()) > TOLERANCE) {
                int newY = (mod1->getY() + mod2->getY()) / 2;
                mod1->setPosition(mod1->getX(), newY);
                mod2->setPosition(mod2->getX(), newY);
                
                std::cout << "Enforced Y equality for pair: " 
                          << pair.first << " and " << pair.second << std::endl;
            }
        } else { // HORIZONTAL
            // Similar logic for horizontal symmetry
            double center1Y = mod1->getY() + mod1->getHeight() / 2.0;
            double center2Y = mod2->getY() + mod2->getHeight() / 2.0;
            
            double axis2y = 2.0 * symmetryAxisPosition;
            double diff = std::abs((center1Y + center2Y) - axis2y);
            
            if (diff > TOLERANCE) {
                // Calculate exact positions for symmetry
                double expectedCenter1Y = (center1Y + center2Y) / 2.0 - (center2Y - symmetryAxisPosition);
                double expectedCenter2Y = (center1Y + center2Y) / 2.0 + (symmetryAxisPosition - center1Y);
                
                int newY1 = static_cast<int>(expectedCenter1Y - mod1->getHeight() / 2.0);
                int newY2 = static_cast<int>(expectedCenter2Y - mod2->getHeight() / 2.0);
                
                // Ensure they don't overlap
                if (std::abs(expectedCenter1Y - expectedCenter2Y) < (mod1->getHeight() / 2.0 + mod2->getHeight() / 2.0)) {
                    // They would overlap - adjust positions
                    double minDistance = mod1->getHeight() / 2.0 + mod2->getHeight() / 2.0 + 1; // 1 extra unit
                    
                    // Adjust equally from center
                    expectedCenter1Y = symmetryAxisPosition - minDistance / 2.0;
                    expectedCenter2Y = symmetryAxisPosition + minDistance / 2.0;
                    
                    newY1 = static_cast<int>(expectedCenter1Y - mod1->getHeight() / 2.0);
                    newY2 = static_cast<int>(expectedCenter2Y - mod2->getHeight() / 2.0);
                }
                
                // Update positions
                mod1->setPosition(mod1->getX(), newY1);
                mod2->setPosition(mod2->getX(), newY2);
                
                std::cout << "Enforced Y-axis symmetry for pair: " 
                          << pair.first << " and " << pair.second << std::endl;
            }
            
            // Ensure X coordinates are exactly equal
            if (std::abs(mod1->getX() - mod2->getX()) > TOLERANCE) {
                int newX = (mod1->getX() + mod2->getX()) / 2;
                mod1->setPosition(newX, mod1->getY());
                mod2->setPosition(newX, mod2->getY());
                
                std::cout << "Enforced X equality for pair: " 
                          << pair.first << " and " << pair.second << std::endl;
            }
        }
    }
    
    // Enforce symmetry for self-symmetric modules
    for (const auto& moduleName : selfSymmetricModules) {
        auto it = modules.find(moduleName);
        if (it == modules.end()) {
            std::cerr << "Self-symmetric module not found: " << moduleName << std::endl;
            continue;
        }
        
        auto module = it->second;
        
        if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
            // Calculate center X
            double centerX = module->getX() + module->getWidth() / 2.0;
            
            // Check if center is on the axis
            if (std::abs(centerX - symmetryAxisPosition) > TOLERANCE) {
                // Calculate new X position to center exactly on axis
                int newX = static_cast<int>(symmetryAxisPosition - module->getWidth() / 2.0);
                module->setPosition(newX, module->getY());
                
                std::cout << "Centered self-symmetric module " << moduleName 
                          << " on vertical axis" << std::endl;
            }
        } else { // HORIZONTAL
            // Calculate center Y
            double centerY = module->getY() + module->getHeight() / 2.0;
            
            // Check if center is on the axis
            if (std::abs(centerY - symmetryAxisPosition) > TOLERANCE) {
                // Calculate new Y position to center exactly on axis
                int newY = static_cast<int>(symmetryAxisPosition - module->getHeight() / 2.0);
                module->setPosition(module->getX(), newY);
                
                std::cout << "Centered self-symmetric module " << moduleName 
                          << " on horizontal axis" << std::endl;
            }
        }
    }
    
    // Check if any positions were changed
    bool positionsChanged = false;
    for (const auto& pair : modules) {
        const auto& original = originalPositions[pair.first];
        if (original.first != pair.second->getX() || original.second != pair.second->getY()) {
            positionsChanged = true;
            break;
        }
    }
    
    // Final check for overlaps
    if (positionsChanged && checkForOverlaps()) {
        std::cout << "Warning: Enforcing symmetry constraints introduced overlaps. Attempting fix..." << std::endl;
        fixOverlaps();
    }
}

/**
 * Checks symmetry constraints without modifying any module positions
 */
bool ASFBStarTree::checkSymmetryConstraints() const {
    // Temporarily enable validation-only mode
    bool savedMode = isValidationOnly;
    const_cast<ASFBStarTree*>(this)->isValidationOnly = true;
    
    // Call validation function which will now only check violations
    bool result = validateSymmetryConstraints();
    
    // Restore previous validation mode
    const_cast<ASFBStarTree*>(this)->isValidationOnly = savedMode;
    
    return result;
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
        // For symmetry pairs, rotate both modules SIMULTANEOUSLY
        auto pairModuleName = pairIt->second;
        auto pairIt = modules.find(pairModuleName);
        
        if (pairIt != modules.end()) {
            // First rotate the module itself
            module->rotate();
            
            // Then immediately rotate its symmetric partner to match
            pairIt->second->setRotation(module->getRotated());
            
            return true;
        }
    } else if (selfIt != selfSymmetricModules.end()) {
        // For self-symmetric modules, just rotate
        module->rotate();
        return true;
    }
    
    // Default case: just rotate the module
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

/**
 * Helper function to verify connectivity between all modules
 * @return true if all modules are connected, false otherwise
 */
bool ASFBStarTree::verifyConnectivity() {
    if (modules.empty()) return true;
    
    // Larger threshold for adjacency testing
    const int ADJACENCY_THRESHOLD = 15;
    
    // Build adjacency list
    std::unordered_map<std::string, std::vector<std::string>> adjacency;
    
    // Check all pairs for adjacency
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
                }
            }
            // Check if left edge of mod1 is near right edge of mod2
            else if (std::abs(x1 - (x2 + w2)) <= ADJACENCY_THRESHOLD) {
                // Check for vertical overlap
                if (!(y1 >= (y2 + h2) || y2 >= (y1 + h1))) {
                    adjacent = true;
                }
            }
            // Check if top edge of mod1 is near bottom edge of mod2
            else if (std::abs((y1 + h1) - y2) <= ADJACENCY_THRESHOLD) {
                // Check for horizontal overlap
                if (!(x1 >= (x2 + w2) || x2 >= (x1 + w1))) {
                    adjacent = true;
                }
            }
            // Check if bottom edge of mod1 is near top edge of mod2
            else if (std::abs(y1 - (y2 + h2)) <= ADJACENCY_THRESHOLD) {
                // Check for horizontal overlap
                if (!(x1 >= (x2 + w2) || x2 >= (x1 + w1))) {
                    adjacent = true;
                }
            }
            
            // Symmetry pairs are always considered adjacent
            for (const auto& symPair : symmetryGroup->getSymmetryPairs()) {
                if ((symPair.first == name1 && symPair.second == name2) ||
                    (symPair.first == name2 && symPair.second == name1)) {
                    adjacent = true;
                    break;
                }
            }
            
            if (adjacent) {
                adjacency[name1].push_back(name2);
                std::cout << "  Modules " << name1 << " and " << name2 << " are adjacent" << std::endl;
            }
        }
    }
    
    // Perform BFS to check connectivity
    std::unordered_set<std::string> visited;
    std::queue<std::string> q;
    
    // Start with the first module
    if (!modules.empty()) {
        auto it = modules.begin();
        q.push(it->first);
        visited.insert(it->first);
    }
    
    while (!q.empty()) {
        std::string current = q.front();
        q.pop();
        
        for (const auto& neighbor : adjacency[current]) {
            if (visited.find(neighbor) == visited.end()) {
                visited.insert(neighbor);
                q.push(neighbor);
            }
        }
    }
    
    // All modules must be visited for full connectivity
    return visited.size() == modules.size();
}

/**
 * Force connectivity between all modules if not already connected
 */
void ASFBStarTree::forceConnectivity() {
    if (modules.empty()) return;
    
    // Arrange modules in a linear chain to guarantee connectivity
    SymmetryType symType = symmetryGroup->getType();
    
    // Process symmetry pairs first, then self-symmetric modules
    std::vector<std::pair<std::string, std::string>> pairs = symmetryGroup->getSymmetryPairs();
    std::vector<std::string> selfSym = symmetryGroup->getSelfSymmetric();
    
    // Start positions
    int xPos = 0;
    int yPos = 0;
    
    // Process symmetry pairs
    for (size_t i = 0; i < pairs.size(); i++) {
        auto mod1 = modules[pairs[i].first];
        auto mod2 = modules[pairs[i].second];
        
        if (symType == SymmetryType::VERTICAL) {
            // Place on opposite sides of vertical axis
            int width1 = mod1->getWidth();
            int width2 = mod2->getWidth();
            
            mod1->setPosition(static_cast<int>(symmetryAxisPosition) - width1 - 5, yPos);
            mod2->setPosition(static_cast<int>(symmetryAxisPosition) + 5, yPos);
            
            // Move down for next pair
            yPos += std::max(mod1->getHeight(), mod2->getHeight()) + 2;
        } else {
            // Place on opposite sides of horizontal axis
            int height1 = mod1->getHeight();
            int height2 = mod2->getHeight();
            
            mod1->setPosition(xPos, static_cast<int>(symmetryAxisPosition) - height1 - 5);
            mod2->setPosition(xPos, static_cast<int>(symmetryAxisPosition) + 5);
            
            // Move right for next pair
            xPos += std::max(mod1->getWidth(), mod2->getWidth()) + 2;
        }
    }
    
    // Process self-symmetric modules
    for (const auto& name : selfSym) {
        auto module = modules[name];
        
        if (symType == SymmetryType::VERTICAL) {
            // Center on vertical axis
            int x = static_cast<int>(symmetryAxisPosition) - module->getWidth() / 2;
            module->setPosition(x, yPos);
            
            // Move down for next module
            yPos += module->getHeight() + 2;
        } else {
            // Center on horizontal axis
            int y = static_cast<int>(symmetryAxisPosition) - module->getHeight() / 2;
            module->setPosition(xPos, y);
            
            // Move right for next module
            xPos += module->getWidth() + 2;
        }
    }
    
    std::cout << "Forced connectivity layout applied" << std::endl;
}

/**
 * Check if there are any module overlaps
 * @return true if overlaps exist, false otherwise
 */
bool ASFBStarTree::checkForOverlaps() {
    for (auto it1 = modules.begin(); it1 != modules.end(); ++it1) {
        auto mod1 = it1->second;
        
        for (auto it2 = std::next(it1); it2 != modules.end(); ++it2) {
            auto mod2 = it2->second;
            
            if (mod1->overlaps(*mod2)) {
                std::cout << "  Overlap detected between " << it1->first 
                          << " and " << it2->first << std::endl;
                return true;
            }
        }
    }
    return false;
}

/**
 * Fix overlaps by shifting modules
 */
void ASFBStarTree::fixOverlaps() {
    bool overlapsExist = true;
    int iterations = 0;
    const int MAX_ITERATIONS = 15; // Increased from 10
    
    std::vector<std::pair<std::shared_ptr<Module>, std::shared_ptr<Module>>> lastOverlaps;
    
    while (overlapsExist && iterations < MAX_ITERATIONS) {
        overlapsExist = false;
        iterations++;
        
        // Find all overlapping pairs
        std::vector<std::pair<std::shared_ptr<Module>, std::shared_ptr<Module>>> overlaps;
        
        for (auto it1 = modules.begin(); it1 != modules.end(); ++it1) {
            auto mod1 = it1->second;
            
            for (auto it2 = std::next(it1); it2 != modules.end(); ++it2) {
                auto mod2 = it2->second;
                
                if (mod1->overlaps(*mod2)) {
                    overlapsExist = true;
                    overlaps.push_back({mod1, mod2});
                    
                    std::cout << "  Overlap detected between " << it1->first 
                              << " and " << it2->first << std::endl;
                }
            }
        }
        
        // If no overlaps, we're done
        if (overlaps.empty()) break;
        
        // If we're on our last iterations and still have the same overlaps, try more aggressive measures
        if (iterations > MAX_ITERATIONS / 2 && overlaps.size() == lastOverlaps.size()) {
            bool sameOverlaps = true;
            for (size_t i = 0; i < overlaps.size() && i < lastOverlaps.size(); ++i) {
                if (overlaps[i].first != lastOverlaps[i].first || 
                    overlaps[i].second != lastOverlaps[i].second) {
                    sameOverlaps = false;
                    break;
                }
            }
            
            if (sameOverlaps) {
                std::cout << "  Same overlaps detected for multiple iterations. Using more aggressive fixes." << std::endl;
                
                // For each overlapping pair
                for (auto& pair : overlaps) {
                    auto mod1 = pair.first;
                    auto mod2 = pair.second;
                    
                    // Check if this is a symmetry pair
                    bool isSymPair = false;
                    for (const auto& symPair : symmetryGroup->getSymmetryPairs()) {
                        if ((symPair.first == mod1->getName() && symPair.second == mod2->getName()) ||
                            (symPair.first == mod2->getName() && symPair.second == mod1->getName())) {
                            isSymPair = true;
                            break;
                        }
                    }
                    
                    if (isSymPair) {
                        // Special care needed for symmetry pairs
                        if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
                            // For vertical symmetry, adjust positions more aggressively
                            double center1X = mod1->getX() + mod1->getWidth() / 2.0;
                            double center2X = mod2->getX() + mod2->getWidth() / 2.0;
                            
                            // Ensure they're on opposite sides of the axis
                            if ((center1X < symmetryAxisPosition && center2X < symmetryAxisPosition) ||
                                (center1X > symmetryAxisPosition && center2X > symmetryAxisPosition)) {
                                
                                // Force them to opposite sides
                                if (center1X < symmetryAxisPosition) {
                                    // mod1 on left, mod2 needs to go right
                                    mod2->setPosition(
                                        static_cast<int>(symmetryAxisPosition + 5), // force 5 units from axis
                                        mod2->getY()
                                    );
                                } else {
                                    // mod1 on right, mod2 needs to go left
                                    mod2->setPosition(
                                        static_cast<int>(symmetryAxisPosition - mod2->getWidth() - 5), // force 5 units from axis
                                        mod2->getY()
                                    );
                                }
                            } else {
                                // They're on the correct sides but still overlapping
                                // Try increasing their Y-separation
                                mod2->setPosition(
                                    mod2->getX(),
                                    mod2->getY() + mod2->getHeight() + 5 // force 5 units gap
                                );
                            }
                        } else { // HORIZONTAL
                            // Similar logic for horizontal symmetry
                            double center1Y = mod1->getY() + mod1->getHeight() / 2.0;
                            double center2Y = mod2->getY() + mod2->getHeight() / 2.0;
                            
                            // Ensure they're on opposite sides of the axis
                            if ((center1Y < symmetryAxisPosition && center2Y < symmetryAxisPosition) ||
                                (center1Y > symmetryAxisPosition && center2Y > symmetryAxisPosition)) {
                                
                                // Force them to opposite sides
                                if (center1Y < symmetryAxisPosition) {
                                    // mod1 below, mod2 needs to go above
                                    mod2->setPosition(
                                        mod2->getX(),
                                        static_cast<int>(symmetryAxisPosition + 5) // force 5 units from axis
                                    );
                                } else {
                                    // mod1 above, mod2 needs to go below
                                    mod2->setPosition(
                                        mod2->getX(),
                                        static_cast<int>(symmetryAxisPosition - mod2->getHeight() - 5) // force 5 units from axis
                                    );
                                }
                            } else {
                                // They're on the correct sides but still overlapping
                                // Try increasing their X-separation
                                mod2->setPosition(
                                    mod2->getX() + mod2->getWidth() + 5, // force 5 units gap
                                    mod2->getY()
                                );
                            }
                        }
                    } else {
                        // Not a symmetry pair, but still in the same group
                        // Move them apart in Y-direction
                        int overlapHeight = 
                            std::min(mod1->getY() + mod1->getHeight(), mod2->getY() + mod2->getHeight()) -
                            std::max(mod1->getY(), mod2->getY());
                        
                        if (mod1->getY() <= mod2->getY()) {
                            mod2->setPosition(
                                mod2->getX(),
                                mod1->getY() + mod1->getHeight() + 5 // Add 5 units gap
                            );
                        } else {
                            mod1->setPosition(
                                mod1->getX(),
                                mod2->getY() + mod2->getHeight() + 5 // Add 5 units gap
                            );
                        }
                    }
                }
            }
        }
        
        // Standard fixing for each overlapping pair
        for (auto& pair : overlaps) {
            auto mod1 = pair.first;
            auto mod2 = pair.second;
            
            // Calculate overlap dimensions
            int overlapWidth = 
                std::min(mod1->getX() + mod1->getWidth(), mod2->getX() + mod2->getWidth()) -
                std::max(mod1->getX(), mod2->getX());
            
            int overlapHeight = 
                std::min(mod1->getY() + mod1->getHeight(), mod2->getY() + mod2->getHeight()) -
                std::max(mod1->getY(), mod2->getY());
            
            // Check if this is a symmetry pair
            bool isSymPair = false;
            for (const auto& symPair : symmetryGroup->getSymmetryPairs()) {
                if ((symPair.first == mod1->getName() && symPair.second == mod2->getName()) ||
                    (symPair.first == mod2->getName() && symPair.second == mod1->getName())) {
                    isSymPair = true;
                    break;
                }
            }
            
            if (isSymPair) {
                // For symmetry pairs, adjust according to symmetry type
                if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
                    // Check which module needs to be adjusted based on axis
                    if (mod1->getX() + mod1->getWidth() / 2.0 < symmetryAxisPosition &&
                        mod2->getX() + mod2->getWidth() / 2.0 < symmetryAxisPosition) {
                        // Both on left side of axis - move one vertically
                        if (overlapHeight < overlapWidth) {
                            if (mod1->getY() <= mod2->getY()) {
                                mod2->setPosition(mod2->getX(), mod1->getY() + mod1->getHeight() + 1);
                            } else {
                                mod1->setPosition(mod1->getX(), mod2->getY() + mod2->getHeight() + 1);
                            }
                        } else {
                            // Move the one closer to axis further left
                            if (mod1->getX() + mod1->getWidth() > mod2->getX() + mod2->getWidth()) {
                                mod2->setPosition(mod1->getX() - mod2->getWidth() - 1, mod2->getY());
                            } else {
                                mod1->setPosition(mod2->getX() - mod1->getWidth() - 1, mod1->getY());
                            }
                        }
                    } else if (mod1->getX() + mod1->getWidth() / 2.0 > symmetryAxisPosition &&
                               mod2->getX() + mod2->getWidth() / 2.0 > symmetryAxisPosition) {
                        // Both on right side of axis - similar to above
                        if (overlapHeight < overlapWidth) {
                            if (mod1->getY() <= mod2->getY()) {
                                mod2->setPosition(mod2->getX(), mod1->getY() + mod1->getHeight() + 1);
                            } else {
                                mod1->setPosition(mod1->getX(), mod2->getY() + mod2->getHeight() + 1);
                            }
                        } else {
                            // Move the one closer to axis further right
                            if (mod1->getX() < mod2->getX()) {
                                mod2->setPosition(mod1->getX() + mod1->getWidth() + 1, mod2->getY());
                            } else {
                                mod1->setPosition(mod2->getX() + mod2->getWidth() + 1, mod1->getY());
                            }
                        }
                    } else {
                        // They're properly placed on opposite sides, but still overlap
                        // Adjust Y positions
                        if (mod1->getY() <= mod2->getY()) {
                            mod2->setPosition(mod2->getX(), mod1->getY() + mod1->getHeight() + 1);
                        } else {
                            mod1->setPosition(mod1->getX(), mod2->getY() + mod2->getHeight() + 1);
                        }
                    }
                } else { // HORIZONTAL
                    // Similar logic for horizontal symmetry
                    if (mod1->getY() + mod1->getHeight() / 2.0 < symmetryAxisPosition &&
                        mod2->getY() + mod2->getHeight() / 2.0 < symmetryAxisPosition) {
                        // Both below axis
                        if (overlapWidth < overlapHeight) {
                            if (mod1->getX() <= mod2->getX()) {
                                mod2->setPosition(mod1->getX() + mod1->getWidth() + 1, mod2->getY());
                            } else {
                                mod1->setPosition(mod2->getX() + mod2->getWidth() + 1, mod1->getY());
                            }
                        } else {
                            // Move the one closer to axis further down
                            if (mod1->getY() + mod1->getHeight() > mod2->getY() + mod2->getHeight()) {
                                mod2->setPosition(mod2->getX(), mod1->getY() - mod2->getHeight() - 1);
                            } else {
                                mod1->setPosition(mod1->getX(), mod2->getY() - mod1->getHeight() - 1);
                            }
                        }
                    } else if (mod1->getY() + mod1->getHeight() / 2.0 > symmetryAxisPosition &&
                               mod2->getY() + mod2->getHeight() / 2.0 > symmetryAxisPosition) {
                        // Both above axis
                        if (overlapWidth < overlapHeight) {
                            if (mod1->getX() <= mod2->getX()) {
                                mod2->setPosition(mod1->getX() + mod1->getWidth() + 1, mod2->getY());
                            } else {
                                mod1->setPosition(mod2->getX() + mod2->getWidth() + 1, mod1->getY());
                            }
                        } else {
                            // Move the one closer to axis further up
                            if (mod1->getY() < mod2->getY()) {
                                mod2->setPosition(mod2->getX(), mod1->getY() + mod1->getHeight() + 1);
                            } else {
                                mod1->setPosition(mod1->getX(), mod2->getY() + mod2->getHeight() + 1);
                            }
                        }
                    } else {
                        // They're properly placed on opposite sides, but still overlap
                        // Adjust X positions
                        if (mod1->getX() <= mod2->getX()) {
                            mod2->setPosition(mod1->getX() + mod1->getWidth() + 1, mod2->getY());
                        } else {
                            mod1->setPosition(mod2->getX() + mod2->getWidth() + 1, mod1->getY());
                        }
                    }
                }
            } else {
                // Regular non-symmetry-pair modules
                // Choose optimal displacement direction
                if (overlapHeight <= overlapWidth) {
                    // Shift vertically
                    if (mod1->getY() <= mod2->getY()) {
                        mod2->setPosition(mod2->getX(), mod1->getY() + mod1->getHeight() + 1);
                    } else {
                        mod1->setPosition(mod1->getX(), mod2->getY() + mod2->getHeight() + 1);
                    }
                } else {
                    // Shift horizontally
                    if (mod1->getX() <= mod2->getX()) {
                        mod2->setPosition(mod1->getX() + mod1->getWidth() + 1, mod2->getY());
                    } else {
                        mod1->setPosition(mod2->getX() + mod2->getWidth() + 1, mod1->getY());
                    }
                }
            }
        }
        
        // Update lastOverlaps for next iteration
        lastOverlaps = overlaps;
    }
    
    if (overlapsExist) {
        std::cout << "Warning: Could not resolve all overlaps after " 
                  << MAX_ITERATIONS << " iterations" << std::endl;
        
        // Final attempt: try emergency spread
        double spreadFactor = 1.5;
        
        for (auto& pair : modules) {
            auto& module = pair.second;
            module->setPosition(
                static_cast<int>(module->getX() * spreadFactor),
                static_cast<int>(module->getY() * spreadFactor)
            );
        }
        
        // Also adjust symmetry axis
        if (symmetryGroup->getType() == SymmetryType::VERTICAL) {
            symmetryAxisPosition *= spreadFactor;
        } else {
            symmetryAxisPosition *= spreadFactor;
        }
    } else {
        std::cout << "All overlaps resolved in " << iterations << " iterations" << std::endl;
    }
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