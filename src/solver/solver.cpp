// solver.cpp
#include "solver.hpp"
#include <iostream>
#include <ctime>
#include <algorithm>
#include <limits>
#include <queue>
#include <random>
#include <cmath>
#include <set>

PlacementSolver::PlacementSolver()
    : hbTree(nullptr),
      initialTemperature(1000.0),
      finalTemperature(0.1),
      coolingRate(0.98), // Slower cooling rate for better convergence
      iterationsPerTemperature(300), // More iterations per temperature
      noImprovementLimit(3000), // Longer no improvement limit
      probRotate(0.4),
      probMove(0.3),
      probSwap(0.2),
      probChangeRep(0.05),
      probConvertSym(0.05),
      areaWeight(1.0),
      wirelengthWeight(0.0),
      randomSeed(static_cast<unsigned int>(std::time(nullptr))),
      totalArea(0),
      solutionWidth(0),
      solutionHeight(0),
      timeLimit(290) { // 4 minutes 50 seconds time limit
}

PlacementSolver::~PlacementSolver() {
    // Smart pointers handle cleanup
}

void PlacementSolver::loadProblem(const std::map<std::string, std::shared_ptr<Module>>& modules,
                                const std::vector<std::shared_ptr<SymmetryGroup>>& symmetryGroups) {
    // Store all modules for reference
    allModules = modules;
    this->symmetryGroups = symmetryGroups;
    
    // unified build
    hbTree = std::make_shared<HBStarTree>();
    
    for (auto &m : modules)            // all modules, symmetric or not
        hbTree->addModule(m.second);
    
    for (auto &g : symmetryGroups)     // group meta already parsed
        hbTree->addSymmetryGroup(g);
    
    hbTree->constructInitialTree();    // left‑skewed or balanced
}

void PlacementSolver::createInitialSolution() {
    if (!hbTree->pack()) {
        std::cerr << "Initial solution packing failed. Attempting recovery..." << std::endl;
        // Try simple recovery by re-initializing
        hbTree->constructInitialTree();
        hbTree->pack();
    }
    
    // Calculate initial cost
    calculateTotalArea();
}

bool PlacementSolver::packSolution() {
    return hbTree->pack();
}

int PlacementSolver::calculateTotalArea() {
    if (!hbTree->isPacked()) {
        if (!hbTree->pack()) {
            std::cerr << "Warning: Packing failed during area calculation" << std::endl;
            return std::numeric_limits<int>::max();
        }
    }
    
    totalArea = hbTree->getArea();
    
    // Update dimensions for statistics
    auto modules = hbTree->getModules();
    int minX = std::numeric_limits<int>::max();
    int minY = std::numeric_limits<int>::max();
    int maxX = 0;
    int maxY = 0;
    
    for (const auto& pair : modules) {
        const auto& module = pair.second;
        
        minX = std::min(minX, module->getX());
        minY = std::min(minY, module->getY());
        maxX = std::max(maxX, module->getX() + module->getWidth());
        maxY = std::max(maxY, module->getY() + module->getHeight());
    }
    
    solutionWidth = maxX - minX;
    solutionHeight = maxY - minY;
    
    return totalArea;
}

bool PlacementSolver::hasOverlaps() const {
    return hbTree->hasOverlap();
}

bool PlacementSolver::validateSymmetryConstraints() const {
    for (const auto& group : symmetryGroups) {
        auto symGroupNode = hbTree->getSymmetryGroupNode(group->getName());
        if (!symGroupNode) continue;
        
        auto asfTree = symGroupNode->getASFTree();
        if (!asfTree) continue;
        
        if (!asfTree->validateSymmetryConstraints()) {
            std::cerr << "Symmetry constraints violated for group: " << group->getName() << std::endl;
            return false;
        }
    }
    
    return true;
}

void PlacementSolver::setAnnealingParameters(double initialTemp, double finalTemp, double coolRate, 
                                           int iterations, int noImprovementLimit) {
    initialTemperature = initialTemp;
    finalTemperature = finalTemp;
    coolingRate = coolRate;
    iterationsPerTemperature = iterations;
    this->noImprovementLimit = noImprovementLimit;
}

void PlacementSolver::setPerturbationProbabilities(double rotate, double move, double swap, 
                                                 double changeRep, double convertSym) {
    // Check if probabilities sum to 1.0
    double sum = rotate + move + swap + changeRep + convertSym;
    if (std::abs(sum - 1.0) > 1e-6) {
        // Normalize probabilities to sum to 1.0
        if (sum <= 0.0) {
            // Default values if all probabilities are zero or negative
            probRotate = 0.4;
            probMove = 0.3;
            probSwap = 0.2;
            probChangeRep = 0.05;
            probConvertSym = 0.05;
            return;
        }
        
        probRotate = rotate / sum;
        probMove = move / sum;
        probSwap = swap / sum;
        probChangeRep = changeRep / sum;
        probConvertSym = convertSym / sum;
    } else {
        probRotate = rotate;
        probMove = move;
        probSwap = swap;
        probChangeRep = changeRep;
        probConvertSym = convertSym;
    }
}

void PlacementSolver::setCostWeights(double area, double wirelength) {
    areaWeight = area;
    wirelengthWeight = wirelength;
}

void PlacementSolver::setRandomSeed(unsigned int seed) {
    randomSeed = seed;
}

void PlacementSolver::setTimeLimit(int seconds) {
    timeLimit = seconds;
}

bool PlacementSolver::solve() {
    // Create initial solution
    createInitialSolution();
    
    if (!hbTree) {
        std::cerr << "Error: No modules to place." << std::endl;
        return false;
    }
    
    // Calculate initial area
    calculateTotalArea();
    std::cout << "Initial area: " << totalArea << std::endl;
    
    // Setup enhanced SA
    std::cout << "Starting simulated annealing..." << std::endl;
    std::cout << "Initial temperature: " << initialTemperature << std::endl;
    std::cout << "Final temperature: " << finalTemperature << std::endl;
    std::cout << "Cooling rate: " << coolingRate << std::endl;
    std::cout << "Iterations per temperature: " << iterationsPerTemperature << std::endl;
    std::cout << "No improvement limit: " << noImprovementLimit << std::endl;
    
    // Save initial solution as best solution
    std::map<std::string, std::pair<int, int>> bestPositions;
    std::map<std::string, bool> bestRotations;
    int bestArea = totalArea;
    
    // Save positions and rotations of all modules
    for (const auto& pair : allModules) {
        const auto& module = pair.second;
        bestPositions[pair.first] = {module->getX(), module->getY()};
        bestRotations[pair.first] = module->getRotated();
    }
    
    // Initialize random number generator
    std::mt19937 rng(randomSeed);
    std::uniform_real_distribution<double> uniformDist(0.0, 1.0);
    
    // SA parameters
    double temperature = initialTemperature;
    int noImprovementCount = 0;
    int totalIterations = 0;
    int acceptedMoves = 0;
    int rejectedMoves = 0;
    
    // Start timing
    auto startTime = std::chrono::steady_clock::now();
    
    // Main SA loop
    while (temperature > finalTemperature && noImprovementCount < noImprovementLimit) {
        // Check time limit
        auto currentTime = std::chrono::steady_clock::now();
        auto elapsedSeconds = std::chrono::duration_cast<std::chrono::seconds>(
            currentTime - startTime).count();
        
        if (elapsedSeconds >= timeLimit) {
            std::cout << "Time limit reached, stopping SA..." << std::endl;
            break;
        }
        
        bool improved = false;
        int acceptedInPassCount = 0;
        
        // Perform iterations at current temperature
        for (int i = 0; i < iterationsPerTemperature; ++i) {
            totalIterations++;
            
            // Save current state before perturbation
            std::map<std::string, std::pair<int, int>> currentPositions;
            std::map<std::string, bool> currentRotations;
            
            // Store the current state of each module
            for (const auto& pair : allModules) {
                const auto& module = pair.second;
                currentPositions[pair.first] = {module->getX(), module->getY()};
                currentRotations[pair.first] = module->getRotated();
            }
            
            // Perform a random perturbation
            bool perturbSuccess = false;
            std::string perturbedModule;
            std::string perturbedGroup;
            
            // Choose perturbation type
            double randVal = uniformDist(rng);
            
            if (randVal < probRotate) {
                // Rotate a random module
                auto nodes = hbTree->getModuleNodes();
                if (!nodes.empty()) {
                    int idx = uniformDist(rng) * nodes.size();
                    auto it = nodes.begin();
                    std::advance(it, idx);
                    perturbedModule = it->first;
                    perturbSuccess = hbTree->rotateModule(perturbedModule);
                }
            } 
            else if (randVal < probRotate + probMove) {
                // Move a node to a new position
                auto nodes = hbTree->getModuleNodes();
                if (nodes.size() >= 2) {
                    int idx1 = uniformDist(rng) * nodes.size();
                    int idx2;
                    do {
                        idx2 = uniformDist(rng) * nodes.size();
                    } while (idx2 == idx1);
                    
                    auto it1 = nodes.begin();
                    auto it2 = nodes.begin();
                    std::advance(it1, idx1);
                    std::advance(it2, idx2);
                    
                    bool asLeftChild = uniformDist(rng) < 0.5;
                    perturbSuccess = hbTree->moveNode(it1->first, it2->first, asLeftChild);
                }
            } 
            else if (randVal < probRotate + probMove + probSwap) {
                // Swap two nodes
                auto nodes = hbTree->getModuleNodes();
                if (nodes.size() >= 2) {
                    int idx1 = uniformDist(rng) * nodes.size();
                    int idx2;
                    do {
                        idx2 = uniformDist(rng) * nodes.size();
                    } while (idx2 == idx1);
                    
                    auto it1 = nodes.begin();
                    auto it2 = nodes.begin();
                    std::advance(it1, idx1);
                    std::advance(it2, idx2);
                    
                    perturbSuccess = hbTree->swapNodes(it1->first, it2->first);
                }
            } 
            else if (randVal < probRotate + probMove + probSwap + probChangeRep) {
                // Change representative
                auto symGroups = hbTree->getSymmetryGroups();
                if (!symGroups.empty()) {
                    int groupIdx = uniformDist(rng) * symGroups.size();
                    auto group = symGroups[groupIdx];
                    
                    const auto& pairs = group->getSymmetryPairs();
                    if (!pairs.empty()) {
                        int pairIdx = uniformDist(rng) * pairs.size();
                        const auto& pair = pairs[pairIdx];
                        
                        perturbedModule = (uniformDist(rng) < 0.5) ? pair.first : pair.second;
                        perturbedGroup = group->getName();
                        
                        perturbSuccess = hbTree->changeRepresentative(perturbedGroup, perturbedModule);
                    }
                }
            } 
            else {
                // Convert symmetry type
                auto symGroups = hbTree->getSymmetryGroups();
                if (!symGroups.empty()) {
                    int groupIdx = uniformDist(rng) * symGroups.size();
                    auto group = symGroups[groupIdx];
                    
                    perturbedGroup = group->getName();
                    perturbSuccess = hbTree->convertSymmetryType(perturbedGroup);
                }
            }
            
            // Skip if perturbation failed
            if (!perturbSuccess) {
                continue;
            }
            
            // Pack solution and check if it's valid
            bool packSuccess = hbTree->pack();
            
            // Check for overlaps and symmetry constraint violations
            bool isValid = packSuccess && !hasOverlaps() && validateSymmetryConstraints();
            
            if (isValid) {
                // Calculate new area
                int newArea = calculateTotalArea();
                int deltaCost = newArea - totalArea;
                
                // Decide whether to accept
                bool accept = false;
                if (deltaCost <= 0) {
                    // Always accept improvements
                    accept = true;
                    if (newArea < bestArea) {
                        improved = true;
                        bestArea = newArea;
                        
                        // Update best solution
                        for (const auto& pair : allModules) {
                            const auto& module = pair.second;
                            bestPositions[pair.first] = {module->getX(), module->getY()};
                            bestRotations[pair.first] = module->getRotated();
                        }
                        
                        // Reset no improvement counter
                        noImprovementCount = 0;
                    }
                } else {
                    // Accept worsening moves with probability based on temperature
                    double acceptProb = std::exp(-deltaCost / temperature);
                    accept = uniformDist(rng) < acceptProb;
                }
                
                if (accept) {
                    // Move accepted
                    totalArea = newArea;
                    acceptedMoves++;
                    acceptedInPassCount++;
                } else {
                    // Move rejected, restore previous state
                    rejectedMoves++;
                    
                    // Restore all module positions and rotations
                    for (auto& pair : allModules) {
                        auto module = pair.second;
                        const auto& oldPos = currentPositions[pair.first];
                        bool oldRot = currentRotations[pair.first];
                        
                        module->setPosition(oldPos.first, oldPos.second);
                        module->setRotation(oldRot);
                    }
                    
                    // Re-pack to ensure internal state is consistent
                    hbTree->pack();
                    
                    // Restore area
                    calculateTotalArea();
                }
            } else {
                // Invalid placement, restore previous state
                rejectedMoves++;
                
                // Restore all module positions and rotations
                for (auto& pair : allModules) {
                    auto module = pair.second;
                    const auto& oldPos = currentPositions[pair.first];
                    bool oldRot = currentRotations[pair.first];
                    
                    module->setPosition(oldPos.first, oldPos.second);
                    module->setRotation(oldRot);
                }
                
                // Re-pack to ensure internal state is consistent
                hbTree->pack();
                
                // Restore area
                calculateTotalArea();
            }
            
            if (!improved) {
                noImprovementCount++;
            }
            
            // Check if we've reached the limit
            if (noImprovementCount >= noImprovementLimit) {
                break;
            }
        }
        
        // Calculate acceptance ratio
        double acceptanceRatio = static_cast<double>(acceptedInPassCount) / iterationsPerTemperature;
        
        // Adaptive cooling schedule
        if (acceptanceRatio > 0.8) {
            // Cool faster if acceptance rate is high
            temperature *= coolingRate * 0.9;
        } else if (acceptanceRatio < 0.1) {
            // Cool slower if acceptance rate is low
            temperature *= coolingRate * 1.1;
            if (temperature > initialTemperature) {
                temperature = initialTemperature;
            }
        } else {
            // Normal cooling
            temperature *= coolingRate;
        }
        
        std::cout << "Temperature: " << temperature 
                  << ", Best area: " << bestArea 
                  << ", Current area: " << totalArea 
                  << ", No improvement: " << noImprovementCount 
                  << ", Acceptance ratio: " << acceptanceRatio 
                  << std::endl;
    }
    
    // Restore best solution
    for (auto& pair : allModules) {
        auto module = pair.second;
        const auto& bestPos = bestPositions[pair.first];
        bool bestRot = bestRotations[pair.first];
        
        module->setPosition(bestPos.first, bestPos.second);
        module->setRotation(bestRot);
    }
    
    // Final pack and validation
    bool finalPackSuccess = hbTree->pack();
    if (!finalPackSuccess) {
        std::cerr << "Final packing failed, trying to recover..." << std::endl;
        
        // Reapply best positions directly
        for (auto& pair : allModules) {
            auto module = pair.second;
            const auto& bestPos = bestPositions[pair.first];
            bool bestRot = bestRotations[pair.first];
            
            module->setPosition(bestPos.first, bestPos.second);
            module->setRotation(bestRot);
        }
        
        // Try packing one more time
        hbTree->pack();
    }
    
    // Final validation
    if (hasOverlaps()) {
        std::cerr << "Warning: Final solution has overlaps" << std::endl;
    }
    
    if (!validateSymmetryConstraints()) {
        std::cerr << "Warning: Final solution violates symmetry constraints" << std::endl;
    }
    
    // Recalculate area
    calculateTotalArea();
    
    std::cout << "Simulated annealing completed." << std::endl;
    std::cout << "Final area: " << totalArea << std::endl;
    
    return true;
}

int PlacementSolver::getSolutionArea() const {
    return totalArea;
}

std::map<std::string, std::shared_ptr<Module>> PlacementSolver::getSolutionModules() const {
    return allModules;
}

std::map<std::string, int> PlacementSolver::getStatistics() const {
    std::map<std::string, int> stats;
    stats["totalArea"] = totalArea;
    stats["width"] = solutionWidth;
    stats["height"] = solutionHeight;
    return stats;
}