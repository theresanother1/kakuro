#ifndef KAKUROGENERATOR_H
#define KAKUROGENERATOR_H

#include <vector>
#include <random>
#include <algorithm>
#include <chrono>
#include "KakuroSolver.h"

// Forward declarations
class KakuroBoard;
class KakuroGenerator;

class KakuroBoard {
private:
    int rows;
    int cols;
    const int minRunLength = 2;
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<> dis;
    std::uniform_int_distribution<> numDis;

    void placeWhiteCells();
    bool validateAndFixRuns();
    bool validateBoard();
    void placeClues();
    std::pair<int, int> calculateRightRun(int row, int col);
    std::pair<int, int> calculateDownRun(int row, int col);

public:
    std::vector<std::vector<Cell>> board;

    KakuroBoard(int r, int c);
    bool generateBoard();
    void printBoard() const;
};

class KakuroGenerator {
private:
    struct Individual {
        std::vector<std::vector<Cell>> board;
        double fitness;

        Individual(int size);
    };

    struct Run {
        int startX, startY;
        int length;
        bool horizontal;

        Run(int x, int y, int l, bool h);
    };

    const int size;
    const int populationSize;
    const double mutationRate;
    const double targetFitness;
    std::mt19937 rng;

    // Constants for fitness calculation
    const double WEIGHT_SOLVABILITY = 0.4;
    const double WEIGHT_UNIQUENESS = 0.3;
    const double WEIGHT_COMPLEXITY = 0.2;
    const double WEIGHT_BALANCE = 0.1;

    double checkBasicValidity(const std::vector<std::vector<Cell>>& board);
    double calculateFitness(Individual& individual);
    Individual createRandomIndividual();
    void heavyMutation(Individual& ind);
    Individual crossover(const Individual& parent1, const Individual& parent2);
    void mutate(Individual& ind, int& generationsWithoutImprovement);
    void printBoard(std::vector<std::vector<Cell>>& board) const;
    bool isClueValuePossible(int clue, int length);
    int countRunLength(const std::vector<std::vector<Cell>>& board, int startX, int startY, bool horizontal);
    void adjustRunOrClue(std::vector<std::vector<Cell>>& board, int x, int y, bool horizontal, std::mt19937& rng);
    void constraintGuidedMutation(Individual& ind, int generationsStuck);
    std::vector<Run> analyzeRunLengths(const std::vector<std::vector<Cell>>& board);
    void optimizeRunLengths(std::vector<std::vector<Cell>>& board);

public:
    KakuroGenerator(int boardSize);
    std::vector<std::vector<Cell>> generateBoard(int maxGenerations = 2000, int maxTimeSeconds = 300);
};

#endif // KAKUROGENERATOR_H