#pragma once

#include <memory>
#include <random>
#include <cmath>
#include <functional>
#include <vector>
#include <string>
#include <chrono>
#include "../data_struct/HBStarTree.hpp"

class SimulatedAnnealing {
private:
    // The current state of the placement solution
    std::shared_ptr<HBStarTree> currentSolution;
    
    // The best solution found so far
    std::shared_ptr<HBStarTree> bestSolution;
    
    // Cost of the current solution
    int currentCost;
    
    // Cost of the best solution
    int bestCost;
    
    vector<string> moduleNameCache;
    vector<string> symGroupNameCache;
    bool cacheInitialized;

    // SA parameters
    double initialTemperature;
    double finalTemperature;
    double coolingRate;
    int iterationsPerTemperature;
    int noImprovementLimit;
    int timeLimit; // Time limit in seconds
    
    // Random number generation - mutable to allow usage in const methods
    mutable std::mt19937 rng;
    mutable std::uniform_real_distribution<double> uniformDist;
    
    // Perturbation probabilities
    double probRotate;
    double probMove;
    double probSwap;
    double probChangeRepresentative;
    double probConvertSymmetryType;
    
    // Statistics
    int totalIterations;
    int acceptedMoves;
    int rejectedMoves;
    int noImprovementCount;
    std::chrono::steady_clock::time_point startTime;
    
    // Cost function weight parameters
    double areaWeight;
    double wirelengthWeight;
    
    // Cache of last perturbation to support incremental updates
    enum class PerturbationType {
        NONE,
        ROTATE,
        MOVE,
        SWAP,
        CHANGE_REP,
        CONVERT_SYM
    };
    
    PerturbationType lastPerturbation;
    std::string lastModule1;
    std::string lastModule2;
    std::string lastSymGroup;
    
    /**
     * Calculates the cost of a solution
     */
    int calculateCost(const std::shared_ptr<HBStarTree>& solution) const;
    
    /**
     * Performs a random perturbation on the current solution
     */
    bool perturb();
    
    /**
     * Performs incremental cost update after perturbation
     * Returns the difference in cost (new - old)
     */
    int updateCostIncremental();
    
    /**
     * Rotates a random module
     */
    bool perturbRotate();
    
    /**
     * Moves a random node to a new position
     */
    bool perturbMove();
    
    /**
     * Swaps two random nodes
     */
    bool perturbSwap();
    
    /**
     * Changes the representative of a symmetry pair in a random symmetry group
     */
    bool perturbChangeRepresentative();
    
    /**
     * Converts the symmetry type of a random symmetry group
     */
    bool perturbConvertSymmetryType();
    
    /**
     * Decides whether to accept a move based on the cost difference and temperature
     * 
     * @param costDifference Difference in cost (new - old)
     * @param temperature Current temperature
     * @return True if the move should be accepted, false otherwise
     */
    bool acceptMove(int costDifference, double temperature) const;
    
    /**
     * Select a random module
     * 
     * @return Name of the selected module
     */
    std::string selectRandomModule() const;
    
    /**
     * Select a random symmetry group
     * 
     * @return Name of the selected symmetry group
     */
    std::string selectRandomSymmetryGroup() const;
    
    /**
     * Select a random node (module or symmetry group)
     * 
     * @return Name of the selected node
     */
    std::string selectRandomNode() const;
    
    /**
     * Checks if time limit is reached
     */
    bool isTimeLimitReached() const;
    
    /**
     * Pre-compute solution information to avoid redundant calculations
     */
    void precomputeSolutionInfo();

public:
    /**
     * Constructor
     */
    SimulatedAnnealing(std::shared_ptr<HBStarTree> initialSolution,
                      double initialTemp = 1000.0,
                      double finalTemp = 0.1,
                      double coolRate = 0.95,
                      int iterations = 100,
                      int noImprovementLimit = 1000,
                      int timeLimit = 300); // 5 minute default

    void initializeNameCache();
    
    /**
     * Sets the perturbation probabilities
     * 
     * @param rotate Probability of rotation operation
     * @param move Probability of move operation
     * @param swap Probability of swap operation
     * @param changeRep Probability of change representative operation
     * @param convertSym Probability of convert symmetry type operation
     */
    void setPerturbationProbabilities(double rotate, double move, double swap, 
                                     double changeRep, double convertSym);
    
    /**
     * Sets the cost function weights
     * 
     * @param area Weight for area term
     * @param wirelength Weight for wirelength term
     */
    void setCostWeights(double area, double wirelength);
    
    /**
     * Runs the simulated annealing algorithm
     * 
     * @return Best solution found
     */
    std::shared_ptr<HBStarTree> run();
    
    std::shared_ptr<HBStarTree> getBestSolution() const;
    int getBestCost() const;
    
    /**
     * Gets statistics about the annealing process
     * 
     * @return Map of statistic name to value
     */
    std::map<std::string, int> getStatistics() const;
    
    /**
     * Sets the random seed for reproducible results
     * 
     * @param seed Random seed
     */
    void setSeed(unsigned int seed);
    
    /**
     * Sets the time limit in seconds
     * 
     * @param seconds Time limit in seconds
     */
    void setTimeLimit(int seconds);
};