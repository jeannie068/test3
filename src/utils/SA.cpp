// SA.cpp
#include "SA.hpp"
#include <iostream>
#include <limits>
#include <queue>
#include <algorithm>
#include <cmath>

SimulatedAnnealing::SimulatedAnnealing(
    const std::map<std::string, std::shared_ptr<Module>>& allModules,
    const std::map<std::string, std::shared_ptr<Module>>& regularModules,
    const std::map<std::string, std::shared_ptr<ASFBStarTree>>& symmetryTrees,
    std::shared_ptr<BStarTreeNode> regularTree,
    const std::map<std::string, std::shared_ptr<SymmetryGroup>>& moduleToGroup,
    double initialTemp,
    double finalTemp,
    double coolRate,
    int iterations,
    int noImprovementLimit,
    int timeLimit
)
    : allModules(allModules),
      regularModules(regularModules),
      symmetryTrees(symmetryTrees),
      regularTree(regularTree),
      moduleToGroup(moduleToGroup),
      initialTemperature(initialTemp),
      finalTemperature(finalTemp),
      coolingRate(coolRate),
      iterationsPerTemperature(iterations),
      noImprovementLimit(noImprovementLimit),
      timeLimit(timeLimit),
      uniformDist(0.0, 1.0),
      // Fix initialization order to match class declaration
      currentArea(0),
      bestArea(std::numeric_limits<int>::max()),
      lastPerturbation(PerturbationType::NONE),
      totalIterations(0),
      acceptedMoves(0),
      rejectedMoves(0),
      noImprovementCount(0),
      probRotate(0.3),
      probMove(0.3),
      probSwap(0.3),
      probChangeRep(0.05),
      probConvertSym(0.05),
      areaWeight(1.0),
      wirelengthWeight(0.0)
{
    // Initialize random number generator
    rng.seed(static_cast<unsigned int>(std::time(nullptr)));
    
    // Pack the initial solution
    if (packSolution()) {
        // Calculate initial cost
        currentArea = calculateArea();
        bestArea = currentArea;
        
        // Save initial solution as best
        saveBestSolution();
    }
}

void SimulatedAnnealing::setPerturbationProbabilities(double rotate, double move, double swap,
                                                    double changeRep, double convertSym) {
    // Normalize probabilities to sum to 1.0
    double sum = rotate + move + swap + changeRep + convertSym;
    if (sum <= 0.0) {
        // Default values
        probRotate = 0.3;
        probMove = 0.3;
        probSwap = 0.3;
        probChangeRep = 0.05;
        probConvertSym = 0.05;
        return;
    }
    
    probRotate = rotate / sum;
    probMove = move / sum;
    probSwap = swap / sum;
    probChangeRep = changeRep / sum;
    probConvertSym = convertSym / sum;
}

void SimulatedAnnealing::setCostWeights(double area, double wirelength) {
    areaWeight = area;
    wirelengthWeight = wirelength;
}

void SimulatedAnnealing::setSeed(unsigned int seed) {
    rng.seed(seed);
}

void SimulatedAnnealing::setTimeLimit(int seconds) {
    timeLimit = seconds;
}

bool SimulatedAnnealing::packSolution() {
    // Start by packing all symmetry groups
    int maxX = 0;
    int maxY = 0;
    
    // First, pack all symmetry islands
    for (auto& pair : symmetryTrees) {
        // Pack the symmetry island
        if (!pair.second->pack()) {
            return false;
        }
        
        // Find the bounding rectangle
        int minX = std::numeric_limits<int>::max();
        int minY = std::numeric_limits<int>::max();
        int symMaxX = 0;
        int symMaxY = 0;
        
        for (const auto& modulePair : pair.second->getModules()) {
            const auto& module = modulePair.second;
            
            minX = std::min(minX, module->getX());
            minY = std::min(minY, module->getY());
            symMaxX = std::max(symMaxX, module->getX() + module->getWidth());
            symMaxY = std::max(symMaxY, module->getY() + module->getHeight());
        }
        
        // Add spacing between symmetry groups to prevent overlaps
        if (maxY > 0) {
            maxY += 5; // Add some spacing between groups
        }
        
        // Shift symmetry island to position (0, maxY)
        int deltaX = -minX;
        int deltaY = maxY - minY;
        
        // Update symmetry axis position
        if (pair.second->getSymmetryGroup()->getType() == SymmetryType::VERTICAL) {
            double oldAxis = pair.second->getSymmetryAxisPosition();
            pair.second->getSymmetryGroup()->setAxisPosition(oldAxis + deltaX);
        } else { // HORIZONTAL
            double oldAxis = pair.second->getSymmetryAxisPosition();
            pair.second->getSymmetryGroup()->setAxisPosition(oldAxis + deltaY);
        }
        
        // Shift all modules in the symmetry island
        for (const auto& modulePair : pair.second->getModules()) {
            auto module = modulePair.second;
            module->setPosition(module->getX() + deltaX, module->getY() + deltaY);
        }
        
        // Update max coordinates
        maxX = std::max(maxX, symMaxX - minX);
        maxY += (symMaxY - minY);
    }
    
    // Then, pack regular modules (if any)
    if (regularTree && !regularModules.empty()) {
        // Add some spacing before regular modules
        if (maxY > 0) {
            maxY += 10;
        }
        
        // Place regular modules
        if (!packRegularModules(0, maxY)) {
            return false;
        }
    }
    
    return !hasOverlaps() && validateSymmetryConstraints();
}

bool SimulatedAnnealing::packRegularModules(int startX, int startY) {
    if (!regularTree) return true; // No regular modules to pack
    
    // Use a simple packing algorithm for B*-tree
    std::queue<std::shared_ptr<BStarTreeNode>> nodeQueue;
    nodeQueue.push(regularTree);
    
    while (!nodeQueue.empty()) {
        auto currentNode = nodeQueue.front();
        nodeQueue.pop();
        
        const std::string& moduleName = currentNode->getModuleName();
        auto moduleIt = regularModules.find(moduleName);
        if (moduleIt == regularModules.end()) {
            continue;
        }
        
        auto module = moduleIt->second;
        
        int x = startX, y = startY;
        
        // Calculate position based on B*-tree rules
        if (currentNode->getParent()) {
            auto parentName = currentNode->getParent()->getModuleName();
            auto parentIt = regularModules.find(parentName);
            if (parentIt != regularModules.end()) {
                auto parentModule = parentIt->second;
                
                if (currentNode->isLeftChild()) {
                    // Left child: place to the right of parent
                    x = parentModule->getX() + parentModule->getWidth();
                    y = parentModule->getY();
                } else {
                    // Right child: place above parent
                    x = parentModule->getX();
                    y = parentModule->getY() + parentModule->getHeight();
                }
            }
        }
        
        // Set the module's position
        module->setPosition(x, y);
        
        // Add children to the queue
        if (currentNode->getLeftChild()) {
            nodeQueue.push(currentNode->getLeftChild());
        }
        if (currentNode->getRightChild()) {
            nodeQueue.push(currentNode->getRightChild());
        }
    }
    
    return true;
}

int SimulatedAnnealing::calculateArea() {
    // Make sure solution is packed before calculating area
    if (!packSolution()) {
        std::cerr << "Warning: Packing failed during area calculation" << std::endl;
        return std::numeric_limits<int>::max();
    }
    
    // Find the bounding rectangle of all modules
    int minX = std::numeric_limits<int>::max();
    int minY = std::numeric_limits<int>::max();
    int maxX = 0;
    int maxY = 0;
    
    // Check regular modules
    for (const auto& pair : allModules) {
        const auto& module = pair.second;
        
        minX = std::min(minX, module->getX());
        minY = std::min(minY, module->getY());
        maxX = std::max(maxX, module->getX() + module->getWidth());
        maxY = std::max(maxY, module->getY() + module->getHeight());
    }
    
    // Also check modules in symmetry trees
    for (const auto& pair : symmetryTrees) {
        for (const auto& modulePair : pair.second->getModules()) {
            const auto& module = modulePair.second;
            
            minX = std::min(minX, module->getX());
            minY = std::min(minY, module->getY());
            maxX = std::max(maxX, module->getX() + module->getWidth());
            maxY = std::max(maxY, module->getY() + module->getHeight());
        }
    }
    
    // Check for invalid dimensions
    if (maxX <= minX || maxY <= minY) {
        std::cerr << "Warning: Invalid dimensions detected in area calculation" << std::endl;
        return std::numeric_limits<int>::max();
    }
    
    // Check for unreasonably large dimensions
    if (maxX > 10000000 || maxY > 10000000) {
        std::cerr << "Warning: Unreasonably large dimensions detected" << std::endl;
        return std::numeric_limits<int>::max();
    }
    
    // Calculate area
    return (maxX - minX) * (maxY - minY);
}


bool SimulatedAnnealing::hasOverlaps() const {
    // Check for overlaps between all pairs of modules
    std::vector<std::shared_ptr<Module>> moduleList;
    
    // Collect all modules
    for (const auto& pair : allModules) {
        moduleList.push_back(pair.second);
    }
    
    // Check all pairs
    for (size_t i = 0; i < moduleList.size(); ++i) {
        for (size_t j = i + 1; j < moduleList.size(); ++j) {
            if (moduleList[i]->overlaps(*moduleList[j])) {
                std::cerr << "Overlap detected between modules: " 
                          << moduleList[i]->getName() << " and " 
                          << moduleList[j]->getName() << std::endl;
                return true;
            }
        }
    }
    
    return false;
}

bool SimulatedAnnealing::validateSymmetryConstraints() const {
    // Check symmetry constraints for all symmetry groups
    for (const auto& pair : symmetryTrees) {
        auto asfTree = pair.second;
        if (!asfTree) continue;
        
        // Use checkSymmetryConstraints (not enforceSymmetryConstraints) to avoid modifying anything
        bool isValid = true;
        
        // First, ensure we're not modifying but just checking
        bool savedMode = asfTree->isValidationOnly;
        const_cast<ASFBStarTree*>(asfTree.get())->isValidationOnly = true;
        
        isValid = asfTree->validateSymmetryConstraints();
        
        // Restore original mode
        const_cast<ASFBStarTree*>(asfTree.get())->isValidationOnly = savedMode;
        
        if (!isValid) {
            std::cerr << "Symmetry constraints violated for group: " << pair.first << std::endl;
            return false;
        }
    }
    
    return true;
}

bool SimulatedAnnealing::perturb() {
    // Reset perturbation tracking variables
    lastPerturbation = PerturbationType::NONE;
    perturbedModule1.clear();
    perturbedModule2.clear();
    perturbedSymGroup.clear();
    
    // Save current state before perturbation
    saveCurrentState();
    
    // Choose perturbation type
    double randVal = uniformDist(rng);
    bool success = false;
    
    if (randVal < probRotate) {
        success = perturbRotate();
        if (success) lastPerturbation = PerturbationType::ROTATE;
    } 
    else if (randVal < probRotate + probMove) {
        success = perturbMove();
        if (success) lastPerturbation = PerturbationType::MOVE;
    } 
    else if (randVal < probRotate + probMove + probSwap) {
        success = perturbSwap();
        if (success) lastPerturbation = PerturbationType::SWAP;
    } 
    else if (randVal < probRotate + probMove + probSwap + probChangeRep) {
        success = perturbChangeRepresentative();
        if (success) lastPerturbation = PerturbationType::CHANGE_REP;
    } 
    else {
        success = perturbConvertSymmetryType();
        if (success) lastPerturbation = PerturbationType::CONVERT_SYM;
    }
    
    return success;
}

bool SimulatedAnnealing::perturbRotate() {
    if (regularModules.empty() && symmetryTrees.empty()) {
        return false;
    }
    
    // Decide whether to rotate a regular module or a module in a symmetry group
    if (!regularModules.empty() && (symmetryTrees.empty() || uniformDist(rng) < 0.5)) {
        // Rotate a regular module
        int index = static_cast<int>(uniformDist(rng) * regularModules.size());
        auto it = regularModules.begin();
        std::advance(it, index);
        perturbedModule1 = it->first;
        
        // Rotate the module
        it->second->rotate();
        return true;
    } 
    else if (!symmetryTrees.empty()) {
        // Rotate a module in a symmetry group
        int index = static_cast<int>(uniformDist(rng) * symmetryTrees.size());
        auto it = symmetryTrees.begin();
        std::advance(it, index);
        perturbedSymGroup = it->first;
        
        // Get a random module from the symmetry group
        auto& modules = it->second->getModules();
        if (modules.empty()) {
            return false;
        }
        
        int moduleIndex = static_cast<int>(uniformDist(rng) * modules.size());
        auto moduleIt = modules.begin();
        std::advance(moduleIt, moduleIndex);
        perturbedModule1 = moduleIt->first;
        
        // Check if this is in a symmetry pair
        auto symGroup = it->second->getSymmetryGroup();
        bool isInSymPair = false;
        for (const auto& pair : symGroup->getSymmetryPairs()) {
            if (pair.first == perturbedModule1 || pair.second == perturbedModule1) {
                isInSymPair = true;
                break;
            }
        }
        
        // For modules in symmetry pairs, we need to rotate both modules
        if (isInSymPair) {
            return it->second->rotateModule(perturbedModule1);
        } else {
            // For self-symmetric modules, just rotate directly
            auto modIt = modules.find(perturbedModule1);
            if (modIt != modules.end()) {
                modIt->second->rotate();
                return true;
            }
        }
    }
    
    return false;
}

bool SimulatedAnnealing::perturbMove() {
    // Implement a simple move perturbation
    // This is a simplification - in a real implementation, you would need to
    // handle moving modules within the B*-tree representation
    
    // For now, let's just try moving a regular module to a new random position
    if (regularModules.empty()) {
        return false;
    }
    
    // Select a random module to move
    int index = static_cast<int>(uniformDist(rng) * regularModules.size());
    auto it = regularModules.begin();
    std::advance(it, index);
    perturbedModule1 = it->first;
    
    // Get current position
    auto module = it->second;
    int oldX = module->getX();
    int oldY = module->getY();
    
    // Generate random offset for movement
    int maxOffset = 50; // Adjust based on your layout size
    int offsetX = static_cast<int>(uniformDist(rng) * maxOffset * 2) - maxOffset;
    int offsetY = static_cast<int>(uniformDist(rng) * maxOffset * 2) - maxOffset;
    
    // Move the module
    module->setPosition(std::max(0, oldX + offsetX), std::max(0, oldY + offsetY));
    
    return true;
}

bool SimulatedAnnealing::perturbSwap() {
    // Implement a swap perturbation
    if (regularModules.size() < 2) {
        return false;
    }
    
    // Select two random modules to swap
    int index1 = static_cast<int>(uniformDist(rng) * regularModules.size());
    int index2;
    do {
        index2 = static_cast<int>(uniformDist(rng) * regularModules.size());
    } while (index1 == index2);
    
    auto it1 = regularModules.begin();
    auto it2 = regularModules.begin();
    std::advance(it1, index1);
    std::advance(it2, index2);
    
    perturbedModule1 = it1->first;
    perturbedModule2 = it2->first;
    
    // Swap positions
    auto module1 = it1->second;
    auto module2 = it2->second;
    
    int tmpX = module1->getX();
    int tmpY = module1->getY();
    
    module1->setPosition(module2->getX(), module2->getY());
    module2->setPosition(tmpX, tmpY);
    
    return true;
}

bool SimulatedAnnealing::perturbChangeRepresentative() {
    // Change representative within a symmetry group
    if (symmetryTrees.empty()) {
        return false;
    }
    
    // Select a random symmetry group
    int groupIndex = static_cast<int>(uniformDist(rng) * symmetryTrees.size());
    auto it = symmetryTrees.begin();
    std::advance(it, groupIndex);
    
    perturbedSymGroup = it->first;
    auto asfTree = it->second;
    auto symGroup = asfTree->getSymmetryGroup();
    
    // Select a random symmetry pair
    if (symGroup->getSymmetryPairs().empty()) {
        return false;
    }
    
    int pairIndex = static_cast<int>(uniformDist(rng) * symGroup->getSymmetryPairs().size());
    auto pair = symGroup->getSymmetryPairs()[pairIndex];
    
    // Randomly choose one module from the pair
    perturbedModule1 = (uniformDist(rng) < 0.5) ? pair.first : pair.second;
    
    // Try to change the representative
    return asfTree->changeRepresentative(perturbedModule1);
}

bool SimulatedAnnealing::perturbConvertSymmetryType() {
    // Convert symmetry type (vertical to horizontal or vice versa)
    if (symmetryTrees.empty()) {
        return false;
    }
    
    // Select a random symmetry group
    int groupIndex = static_cast<int>(uniformDist(rng) * symmetryTrees.size());
    auto it = symmetryTrees.begin();
    std::advance(it, groupIndex);
    
    perturbedSymGroup = it->first;
    auto asfTree = it->second;
    
    // Convert symmetry type
    return asfTree->convertSymmetryType();
}

void SimulatedAnnealing::revertPerturbation() {
    // Restore module positions and rotations
    for (const auto& pair : preperturbationState) {
        auto modIt = allModules.find(pair.first);
        if (modIt != allModules.end()) {
            auto& module = modIt->second;
            const auto& state = pair.second;
            module->setPosition(state.x, state.y);
            module->setRotation(state.rotated);
        }
    }
    
    // Restore symmetry group states
    for (const auto& pair : preperturbationSymState) {
        auto symIt = symmetryTrees.find(pair.first);
        if (symIt != symmetryTrees.end()) {
            auto& tree = symIt->second;
            auto group = tree->getSymmetryGroup();
            const auto& state = pair.second;
            group->setAxisPosition(state.axisPosition);
            group->setType(state.type);
        }
    }
}

bool SimulatedAnnealing::acceptMove(int costDiff, double temp) const {
    // Always accept improvements
    if (costDiff <= 0) {
        return true;
    }
    
    // Use mutable members for random generation in const function
    double probability = std::exp(-costDiff / temp);
    return uniformDist(rng) < probability;
}

void SimulatedAnnealing::saveCurrentState() {
    preperturbationState.clear();
    preperturbationSymState.clear();
    
    // Save state of all modules
    for (const auto& pair : allModules) {
        const auto& module = pair.second;
        preperturbationState[pair.first] = ModuleState(
            module->getX(), module->getY(), module->getRotated());
    }
    
    // Save state of all symmetry groups
    for (const auto& pair : symmetryTrees) {
        const auto& symGroup = pair.second->getSymmetryGroup();
        preperturbationSymState[pair.first] = SymmetryGroupState(
            symGroup->getAxisPosition(), symGroup->getType());
    }
}

void SimulatedAnnealing::saveBestSolution() {
    bestSolution.clear();
    bestSymmetryState.clear();
    
    // Save positions and rotations of all modules
    for (const auto& pair : allModules) {
        const auto& module = pair.second;
        bestSolution[pair.first] = ModuleState(
            module->getX(), module->getY(), module->getRotated());
    }
    
    // Save symmetry group states
    for (const auto& pair : symmetryTrees) {
        const auto& symGroup = pair.second->getSymmetryGroup();
        bestSymmetryState[pair.first] = SymmetryGroupState(
            symGroup->getAxisPosition(), symGroup->getType());
    }
}

void SimulatedAnnealing::restoreBestSolution() {
    // Restore all module positions and rotations
    for (const auto& pair : bestSolution) {
        auto modIt = allModules.find(pair.first);
        if (modIt != allModules.end()) {
            auto& module = modIt->second;
            const auto& state = pair.second;
            
            module->setPosition(state.x, state.y);
            module->setRotation(state.rotated);
        }
    }
    
    // Restore symmetry group states
    for (const auto& pair : bestSymmetryState) {
        auto symIt = symmetryTrees.find(pair.first);
        if (symIt != symmetryTrees.end()) {
            auto& tree = symIt->second;
            auto group = tree->getSymmetryGroup();
            const auto& state = pair.second;
            
            group->setAxisPosition(state.axisPosition);
            group->setType(state.type);
        }
    }
    
    // Update current area to best area
    currentArea = bestArea;
}

bool SimulatedAnnealing::isTimeLimitReached() const {
    auto currentTime = std::chrono::steady_clock::now();
    auto elapsedSeconds = std::chrono::duration_cast<std::chrono::seconds>(
        currentTime - startTime).count();
    return elapsedSeconds >= timeLimit;
}

bool SimulatedAnnealing::run() {
    double temperature = initialTemperature;
    noImprovementCount = 0;
    totalIterations = 0;
    acceptedMoves = 0;
    rejectedMoves = 0;
    startTime = std::chrono::steady_clock::now();
    
    // Pack the initial solution and validate
    if (!packSolution()) {
        std::cerr << "Initial solution is invalid, attempting to fix..." << std::endl;
        
        // Try to fix the initial solution by adjusting positions
        for (auto& pair : symmetryTrees) {
            auto tree = pair.second;
            if (!tree->pack()) {
                std::cerr << "Failed to pack symmetry group: " << pair.first << std::endl;
            }
            
            // IMPORTANT: Enforce symmetry constraints after packing
            tree->enforceSymmetryConstraints();
        }
        
        if (!packSolution() || hasOverlaps()) {
            std::cerr << "Failed to create valid initial solution after fix attempts." << std::endl;
            return false;
        }
    }
    
    // Calculate initial area
    currentArea = calculateArea();
    
    // Save initial solution as best
    bestArea = currentArea;
    saveBestSolution();
    
    std::cout << "Starting SA with initial area: " << currentArea << std::endl;
    
    // Main annealing loop with enhanced cooling and acceptance strategies
    while (temperature > finalTemperature && 
           noImprovementCount < noImprovementLimit && 
           !isTimeLimitReached()) {
        
        int acceptedInPass = 0;
        int totalInPass = 0;
        
        // Iterations at current temperature
        for (int i = 0; i < iterationsPerTemperature && !isTimeLimitReached(); ++i) {
            // Save current state before perturbation
            saveCurrentState();
            
            // Perform perturbation
            if (!perturb()) {
                continue;
            }
            
            totalInPass++;
            
            // Pack the solution
            bool packSuccess = packSolution();
            
            // IMPORTANT CHANGE: Enforce symmetry BEFORE checking validity
            if (packSuccess) {
                // Enforce symmetry constraints for all symmetry groups
                for (auto& pair : symmetryTrees) {
                    pair.second->enforceSymmetryConstraints();
                }
                
                // Now check if the solution is valid with enforced symmetry
                bool isValid = !hasOverlaps() && validateSymmetryConstraints();
                
                if (isValid) {
                    // Calculate new cost
                    int newArea = calculateArea();
                    int costDiff = newArea - currentArea;
                    
                    // Accept or reject with stronger bias towards improvement
                    if (costDiff < 0 || (acceptMove(costDiff, temperature) && newArea < currentArea * 1.2)) {
                        // Accept the perturbation with stricter criteria
                        currentArea = newArea;
                        acceptedMoves++;
                        acceptedInPass++;
                        
                        // Update best solution if improved
                        if (newArea < bestArea) {
                            bestArea = newArea;
                            saveBestSolution();
                            noImprovementCount = 0;
                            
                            std::cout << "New best area found: " << bestArea << std::endl;
                        } else {
                            noImprovementCount++;
                        }
                    } else {
                        // Reject the perturbation
                        revertPerturbation();
                        rejectedMoves++;
                        noImprovementCount++;
                    }
                } else {
                    // Invalid solution, revert the perturbation and try to fix
                    revertPerturbation();
                    
                    // For debugging
                    if (hasOverlaps()) {
                        std::cout << "Placement has overlapping modules" << std::endl;
                    }
                    if (!validateSymmetryConstraints()) {
                        std::cout << "Symmetry constraints violated" << std::endl;
                    }
                    
                    rejectedMoves++;
                }
            } else {
                // Packing failed, revert the perturbation
                revertPerturbation();
                std::cout << "Packing failed" << std::endl;
                rejectedMoves++;
            }
            
            totalIterations++;
            
            // Periodically output stats
            if (totalIterations % 100 == 0) {
                std::cout << "Iteration: " << totalIterations 
                          << ", Temperature: " << temperature 
                          << ", Current area: " << currentArea 
                          << ", Best area: " << bestArea 
                          << std::endl;
            }
        }
        
        // Calculate acceptance ratio with safety check
        double acceptanceRatio = (totalInPass > 0) ? 
            static_cast<double>(acceptedInPass) / totalInPass : 0.0;
        
        // Updated adaptive cooling schedule
        if (acceptanceRatio > 0.6) {
            // Too many acceptances, cool faster
            temperature *= coolingRate * 0.8;
        } else if (acceptanceRatio < 0.1) {
            // Too few acceptances, cool slower
            temperature *= coolingRate * 1.1;
            if (temperature > initialTemperature) {
                temperature = initialTemperature;
            }
        } else {
            // Normal cooling
            temperature *= coolingRate;
        }
        
        // If stuck in a bad local minimum, consider reheating or restarting
        if (acceptanceRatio < 0.05 && noImprovementCount > noImprovementLimit / 2) {
            // Reheat to escape local minimum
            temperature = std::min(initialTemperature * 0.5, temperature * 5.0);
            std::cout << "Reheating to escape local minimum. New temperature: " << temperature << std::endl;
            
            // Consider reverting to best solution to reset path
            if (noImprovementCount > noImprovementLimit * 0.8) {
                restoreBestSolution();
                std::cout << "Restoring best solution to reset path. Area: " << bestArea << std::endl;
            }
        }
        
        std::cout << "Temperature: " << temperature 
                  << ", Best area: " << bestArea 
                  << ", Current area: " << currentArea 
                  << ", No improvement: " << noImprovementCount 
                  << ", Acceptance ratio: " << acceptanceRatio 
                  << std::endl;
    }
    
    // Always restore the best solution at the end
    restoreBestSolution();
    
    // IMPORTANT: Enforce symmetry constraints on the final solution
    for (auto& pair : symmetryTrees) {
        pair.second->enforceSymmetryConstraints();
    }
    
    // Final validation
    bool finalValid = packSolution() && !hasOverlaps() && validateSymmetryConstraints();
    
    if (!finalValid) {
        std::cerr << "Warning: Final solution has issues, attempting to fix..." << std::endl;
        
        // Additional recovery code - try one more time to pack everything
        for (auto& pair : symmetryTrees) {
            auto tree = pair.second;
            if (!tree->pack()) {
                std::cerr << "Failed to pack symmetry group: " << pair.first << std::endl;
            }
            
            // Enforce symmetry constraints
            tree->enforceSymmetryConstraints();
        }
        
        finalValid = packSolution() && !hasOverlaps() && validateSymmetryConstraints();
    }
    
    // If still invalid, try emergency recovery as a last resort
    if (!finalValid) {
        std::cerr << "Warning: Final solution still invalid, attempting emergency recovery..." << std::endl;
        
        for (auto& pair : symmetryTrees) {
            auto tree = pair.second;
            tree->emergencyRecovery();
        }
        
        // Final check after emergency recovery
        finalValid = packSolution() && !hasOverlaps() && validateSymmetryConstraints();
    }
    
    return finalValid;
}

int SimulatedAnnealing::getBestArea() const {
    return bestArea;
}

std::map<std::string, int> SimulatedAnnealing::getStatistics() const {
    std::map<std::string, int> stats;
    stats["totalIterations"] = totalIterations;
    stats["acceptedMoves"] = acceptedMoves;
    stats["rejectedMoves"] = rejectedMoves;
    stats["noImprovementCount"] = noImprovementCount;
    stats["bestArea"] = bestArea;
    
    auto currentTime = std::chrono::steady_clock::now();
    stats["elapsedTimeSeconds"] = static_cast<int>(std::chrono::duration_cast<std::chrono::seconds>(
        currentTime - startTime).count());
    
    return stats;
}