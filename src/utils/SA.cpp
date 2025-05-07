#include "../utils/SA.hpp"
#include <iostream>
#include <algorithm>
#include <ctime>
#include <limits>
#include <cmath>
#include <unordered_set>
#include <chrono>

/**
 * Constructor
 */
SimulatedAnnealing::SimulatedAnnealing(std::shared_ptr<HBStarTree> initialSolution,
                                     double initialTemp,
                                     double finalTemp,
                                     double coolRate,
                                     int iterations,
                                     int noImprovementLimit,
                                     int timeLimit)
    : currentSolution(initialSolution),
      bestSolution(nullptr),
      initialTemperature(initialTemp),
      finalTemperature(finalTemp),
      coolingRate(coolRate),
      iterationsPerTemperature(iterations),
      noImprovementLimit(noImprovementLimit),
      timeLimit(timeLimit),
      uniformDist(0.0, 1.0),
      probRotate(0.3),
      probMove(0.3),
      probSwap(0.3),
      probChangeRepresentative(0.05),
      probConvertSymmetryType(0.05),
      totalIterations(0),
      acceptedMoves(0),
      rejectedMoves(0),
      noImprovementCount(0),
      areaWeight(1.0),
      wirelengthWeight(0.0),
      lastPerturbation(PerturbationType::NONE) {
    
    // Initialize random number generator with current time
    rng.seed(static_cast<unsigned int>(std::time(nullptr)));
    
    // Pack initial solution to get valid coordinates and area
    if (initialSolution && !initialSolution->isPacked()) {
        initialSolution->pack();
    }
    
    // Calculate initial cost
    currentCost = calculateCost(currentSolution);
    
    // Initialize best solution
    bestSolution = currentSolution->clone();
    bestCost = currentCost;
    
    // Precompute solution information
    precomputeSolutionInfo();
}

/**
 * Sets the perturbation probabilities
 */
void SimulatedAnnealing::setPerturbationProbabilities(double rotate, double move, double swap, 
                                                    double changeRep, double convertSym) {
    // Normalize probabilities to sum to 1.0
    double sum = rotate + move + swap + changeRep + convertSym;
    if (sum <= 0.0) {
        // Default values if all probabilities are zero or negative
        probRotate = 0.3;
        probMove = 0.3;
        probSwap = 0.3;
        probChangeRepresentative = 0.05;
        probConvertSymmetryType = 0.05;
        return;
    }
    
    probRotate = rotate / sum;
    probMove = move / sum;
    probSwap = swap / sum;
    probChangeRepresentative = changeRep / sum;
    probConvertSymmetryType = convertSym / sum;
}

/**
 * Sets the cost function weights
 */
void SimulatedAnnealing::setCostWeights(double area, double wirelength) {
    areaWeight = area;
    wirelengthWeight = wirelength;
}

/**
 * Sets the time limit in seconds
 */
void SimulatedAnnealing::setTimeLimit(int seconds) {
    timeLimit = seconds;
}

/**
 * Precompute solution information to avoid redundant calculations
 * This includes things like cache of node connections, etc.
 */
void SimulatedAnnealing::precomputeSolutionInfo() {
    // Here you could precompute and cache information about the solution
    // such as module adjacency lists, etc.
    // This is a placeholder for now
}

/**
 * Calculates the cost of a solution
 */
int SimulatedAnnealing::calculateCost(const std::shared_ptr<HBStarTree>& solution) const {
    if (!solution || !solution->isPacked()) {
        return std::numeric_limits<int>::max();
    }
    
    // Area cost
    int areaCost = solution->getArea();
    
    // Wirelength cost
    int wirelengthCost = solution->getWireLength();
    
    // Weighted sum
    return static_cast<int>(areaWeight * areaCost + wirelengthWeight * wirelengthCost);
}

/**
 * Performs incremental cost update after perturbation
 * Returns the difference in cost (new - old)
 */
int SimulatedAnnealing::updateCostIncremental() {
    // Use incremental packing instead of full packing
    if (!currentSolution->incrementalPack()) {
        // If incremental packing fails, fall back to full packing
        currentSolution->pack();
    }
    
    int newCost = calculateCost(currentSolution);
    int costDifference = newCost - currentCost;
    
    return costDifference;
}

/**
 * Performs a random perturbation on the current solution
 */
bool SimulatedAnnealing::perturb() {
    // Choose perturbation type based on probabilities
    double randVal = uniformDist(rng);
    
    bool success = false;
    
    // Clear last perturbation info
    lastPerturbation = PerturbationType::NONE;
    lastModule1.clear();
    lastModule2.clear();
    lastSymGroup.clear();
    
    if (randVal < probRotate) {
        success = perturbRotate();
        if (success) lastPerturbation = PerturbationType::ROTATE;
    } else if (randVal < probRotate + probMove) {
        success = perturbMove();
        if (success) lastPerturbation = PerturbationType::MOVE;
    } else if (randVal < probRotate + probMove + probSwap) {
        success = perturbSwap();
        if (success) lastPerturbation = PerturbationType::SWAP;
    } else if (randVal < probRotate + probMove + probSwap + probChangeRepresentative) {
        success = perturbChangeRepresentative();
        if (success) lastPerturbation = PerturbationType::CHANGE_REP;
    } else {
        success = perturbConvertSymmetryType();
        if (success) lastPerturbation = PerturbationType::CONVERT_SYM;
    }
    
    return success;
}

/**
 * Rotates a random module
 */
bool SimulatedAnnealing::perturbRotate() {
    std::string moduleName = selectRandomModule();
    if (moduleName.empty()) return false;
    
    lastModule1 = moduleName;
    return currentSolution->rotateModule(moduleName);
}

/**
 * Moves a random node to a new position
 */
bool SimulatedAnnealing::perturbMove() {
    std::string nodeName = selectRandomNode();
    std::string newParentName = selectRandomNode();
    
    if (nodeName.empty() || newParentName.empty() || nodeName == newParentName) {
        return false;
    }
    
    // Randomly decide if the node should be a left or right child
    bool asLeftChild = (uniformDist(rng) < 0.5);
    
    lastModule1 = nodeName;
    lastModule2 = newParentName;
    return currentSolution->moveNode(nodeName, newParentName, asLeftChild);
}

/**
 * Swaps two random nodes
 */
bool SimulatedAnnealing::perturbSwap() {
    std::string nodeName1 = selectRandomNode();
    std::string nodeName2 = selectRandomNode();
    
    if (nodeName1.empty() || nodeName2.empty() || nodeName1 == nodeName2) {
        return false;
    }
    
    lastModule1 = nodeName1;
    lastModule2 = nodeName2;
    return currentSolution->swapNodes(nodeName1, nodeName2);
}

/**
 * Changes the representative of a symmetry pair in a random symmetry group
 */
bool SimulatedAnnealing::perturbChangeRepresentative() {
    std::string symmetryGroupName = selectRandomSymmetryGroup();
    if (symmetryGroupName.empty()) return false;
    
    // Get a random module from the symmetry group
    const auto& symmetryGroups = currentSolution->getSymmetryGroups();
    auto it = std::find_if(symmetryGroups.begin(), symmetryGroups.end(),
                          [&symmetryGroupName](const std::shared_ptr<SymmetryGroup>& group) {
                              return group->getName() == symmetryGroupName;
                          });
    
    if (it == symmetryGroups.end() || (*it)->getSymmetryPairs().empty()) {
        return false;
    }
    
    // Choose a random symmetry pair
    const auto& pairs = (*it)->getSymmetryPairs();
    std::uniform_int_distribution<int> pairDist(0, pairs.size() - 1);
    const auto& pair = pairs[pairDist(rng)];
    
    // Randomly choose one of the modules in the pair
    std::string moduleName = (uniformDist(rng) < 0.5) ? pair.first : pair.second;
    
    lastSymGroup = symmetryGroupName;
    lastModule1 = moduleName;
    return currentSolution->changeRepresentative(symmetryGroupName, moduleName);
}

/**
 * Converts the symmetry type of a random symmetry group
 */
bool SimulatedAnnealing::perturbConvertSymmetryType() {
    std::string symmetryGroupName = selectRandomSymmetryGroup();
    if (symmetryGroupName.empty()) return false;
    
    lastSymGroup = symmetryGroupName;
    return currentSolution->convertSymmetryType(symmetryGroupName);
}

/**
 * Select a random module
 */
std::string SimulatedAnnealing::selectRandomModule() const {
    const auto& modules = currentSolution->getModules();
    if (modules.empty()) return "";
    
    // Convert map to vector for random selection
    std::vector<std::string> moduleNames;
    moduleNames.reserve(modules.size());
    for (const auto& pair : modules) {
        moduleNames.push_back(pair.first);
    }
    
    std::uniform_int_distribution<int> dist(0, moduleNames.size() - 1);
    return moduleNames[dist(rng)];
}

/**
 * Select a random symmetry group
 */
std::string SimulatedAnnealing::selectRandomSymmetryGroup() const {
    const auto& symmetryGroups = currentSolution->getSymmetryGroups();
    if (symmetryGroups.empty()) return "";
    
    std::uniform_int_distribution<int> dist(0, symmetryGroups.size() - 1);
    return symmetryGroups[dist(rng)]->getName();
}

/**
 * Select a random node (module or symmetry group)
 */
std::string SimulatedAnnealing::selectRandomNode() const {
    // Decide whether to select a module or a symmetry group
    const auto& modules = currentSolution->getModules();
    const auto& symmetryGroups = currentSolution->getSymmetryGroups();
    
    int totalNodes = modules.size() + symmetryGroups.size();
    if (totalNodes == 0) return "";
    
    std::uniform_int_distribution<int> dist(0, totalNodes - 1);
    int index = dist(rng);
    
    if (index < static_cast<int>(modules.size())) {
        // Select a module
        auto it = modules.begin();
        std::advance(it, index);
        return it->first;
    } else {
        // Select a symmetry group
        index -= modules.size();
        return symmetryGroups[index]->getName();
    }
}

/**
 * Decides whether to accept a move based on the cost difference and temperature
 */
bool SimulatedAnnealing::acceptMove(int costDifference, double temperature) const {
    // Always accept moves that improve the solution
    if (costDifference <= 0) {
        return true;
    }
    
    // For moves that worsen the solution, accept with a probability based on temperature
    double probability = std::exp(-costDifference / temperature);
    return uniformDist(rng) < probability;
}

/**
 * Checks if time limit is reached
 */
bool SimulatedAnnealing::isTimeLimitReached() const {
    auto currentTime = std::chrono::steady_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::seconds>(
        currentTime - startTime).count();
    return elapsedTime >= timeLimit;
}

/**
 * Runs the SA
 */
std::shared_ptr<HBStarTree> SimulatedAnnealing::run() {
    double temperature = initialTemperature;
    
    // Reset statistics
    totalIterations = 0;
    acceptedMoves = 0;
    rejectedMoves = 0;
    noImprovementCount = 0;
    startTime = std::chrono::steady_clock::now();
    
    // Variable to track the acceptance ratio
    int acceptedInCurrentPass = 0;
    int totalInCurrentPass = 0;
    
    // Main annealing loop
    while (temperature > finalTemperature && 
           noImprovementCount < noImprovementLimit && 
           !isTimeLimitReached()) {
        
        acceptedInCurrentPass = 0;
        totalInCurrentPass = 0;
        
        // Perform iterations at current temperature
        for (int i = 0; i < iterationsPerTemperature && !isTimeLimitReached(); ++i) {
            // Perform a perturbation
            bool perturbSuccess = perturb();
            if (!perturbSuccess) {
                continue;  // Skip this iteration if perturbation failed
            }
            
            totalInCurrentPass++;
            totalIterations++;
            
            // Calculate new cost incrementally
            int costDifference = updateCostIncremental();
            
            // Decide whether to accept the move
            if (acceptMove(costDifference, temperature)) {
                // Accept the move
                currentCost += costDifference;
                acceptedMoves++;
                acceptedInCurrentPass++;
                
                // Update best solution if improved
                if (currentCost < bestCost) {
                    if (!bestSolution) {
                        bestSolution = currentSolution->clone();
                    } else {
                        // Clear previous best solution to avoid memory leak
                        bestSolution.reset();
                        bestSolution = currentSolution->clone();
                    }
                    bestCost = currentCost;
                    noImprovementCount = 0;
                } else {
                    noImprovementCount++;
                }
            } else {
                // Reject the move - revert the perturbation
                switch (lastPerturbation) {
                    case PerturbationType::ROTATE:
                        // Rotate back the module
                        currentSolution->rotateModule(lastModule1);
                        break;
                    case PerturbationType::MOVE:
                        // For move operations, it's complex to revert directly
                        // Just do a full pack to restore consistency
                        currentSolution->pack();
                        break;
                    case PerturbationType::SWAP:
                        // Swap the nodes back
                        currentSolution->swapNodes(lastModule1, lastModule2);
                        break;
                    case PerturbationType::CHANGE_REP:
                        // Change representative back
                        currentSolution->changeRepresentative(lastSymGroup, lastModule1);
                        break;
                    case PerturbationType::CONVERT_SYM:
                        // Convert symmetry type back
                        currentSolution->convertSymmetryType(lastSymGroup);
                        break;
                    case PerturbationType::NONE:
                        // Nothing to revert
                        break;
                }
                
                rejectedMoves++;
                noImprovementCount++;
            }
            
            // Provide progress feedback occasionally
            if (totalIterations % 1000 == 0) {
                std::cout << "Iterations: " << totalIterations 
                          << ", Temp: " << temperature 
                          << ", Best cost: " << bestCost 
                          << ", Current cost: " << currentCost 
                          << std::endl;
            }
        }
        
        // Calculate acceptance ratio
        double acceptanceRatio = (totalInCurrentPass > 0) ? 
            static_cast<double>(acceptedInCurrentPass) / totalInCurrentPass : 0.0;
        
        // Adaptive cooling based on acceptance ratio
        if (acceptanceRatio > 0.8) {
            // Too many acceptances, cool down faster
            temperature *= (coolingRate * 0.9);
        } else if (acceptanceRatio < 0.2) {
            // Too few acceptances, cool down slower
            temperature *= (coolingRate * 1.1);
            if (temperature > initialTemperature) {
                temperature = initialTemperature; // Prevent temperature from exceeding initial
            }
        } else {
            // Normal cooling
            temperature *= coolingRate;
        }
        
        // Print progress
        std::cout << "Temperature: " << temperature 
                  << ", Best cost: " << bestCost 
                  << ", Current cost: " << currentCost 
                  << ", No improvement: " << noImprovementCount 
                  << ", Acceptance ratio: " << acceptanceRatio
                  << std::endl;
                  
        // Early stopping if we have a good solution and low acceptance ratio
        if (acceptanceRatio < 0.01 && noImprovementCount > noImprovementLimit / 2) {
            std::cout << "Early stopping: Low acceptance ratio and no improvement" << std::endl;
            break;
        }
    }
    
    // Make sure the best solution is packed
    if (bestSolution && !bestSolution->isPacked()) {
        bestSolution->pack();
    }
    
    // Report the reason for termination
    if (temperature <= finalTemperature) {
        std::cout << "Annealing finished: Reached final temperature." << std::endl;
    } else if (noImprovementCount >= noImprovementLimit) {
        std::cout << "Annealing finished: Reached no improvement limit." << std::endl;
    } else if (isTimeLimitReached()) {
        std::cout << "Annealing finished: Reached time limit." << std::endl;
    }
    
    return bestSolution;
}

std::shared_ptr<HBStarTree> SimulatedAnnealing::getBestSolution() const {
    return bestSolution;
}

int SimulatedAnnealing::getBestCost() const {
    return bestCost;
}

/**
 * Gets statistics about the annealing process
 */
std::map<std::string, int> SimulatedAnnealing::getStatistics() const {
    std::map<std::string, int> stats;
    stats["totalIterations"] = totalIterations;
    stats["acceptedMoves"] = acceptedMoves;
    stats["rejectedMoves"] = rejectedMoves;
    stats["noImprovementCount"] = noImprovementCount;
    
    // Calculate elapsed time
    auto currentTime = std::chrono::steady_clock::now();
    stats["elapsedTimeSeconds"] = static_cast<int>(std::chrono::duration_cast<std::chrono::seconds>(
        currentTime - startTime).count());
    
    return stats;
}

/**
 * Sets the random seed for reproducible results
 */
void SimulatedAnnealing::setSeed(unsigned int seed) {
    rng.seed(seed);
}