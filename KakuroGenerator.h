//
// Created by Theresa on 19.01.2025.
//

#ifndef PART2_KAKUROGENERATOR_H
#define PART2_KAKUROGENERATOR_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <iomanip>
#include <string>
#include <functional>
#include <chrono>
#include <fstream>
#include <random>
#include "KakuroSolver.h"


class KakuroGenerator {
private:
    int size;
    int populationSize;
    int generations;
    double mutationRate;
    double crossoverRate;
    std::mt19937 rng;

    struct Individual {
        std::vector<std::vector<Cell>> board;
        double fitness;

        Individual(int size) : board(size, std::vector<Cell>(size)), fitness(0.0) {}
    };

    // Helper function to generate a random board layout
    std::vector<std::vector<Cell>> generateRandomBoard();

    // Generate a valid clue for a given run length
    int generateValidClue(int length);

    // Calculate fitness of a board
    double calculateFitness(const std::vector<std::vector<Cell>> &board, int &generation);

    // Evaluate connectivity of white cells
    double evaluateConnectivity(const std::vector<std::vector<Cell>> &board);

    // DFS helper for connectivity check
    void dfs(const std::vector<std::vector<Cell>> &board,
             std::vector<std::vector<bool>> &visited, int i, int j);

    // Evaluate distribution of clues
    double evaluateClueDistribution(const std::vector<std::vector<Cell>> &board);

    // Check if a clue is valid for a given length
    bool isValidClue(int sum, int length);

    // Crossover two boards to create a new one
    std::vector<std::vector<Cell>> crossover(const std::vector<std::vector<Cell>> &parent1,
                                             const std::vector<std::vector<Cell>> &parent2);

    // Mutate a board
    void mutate(std::vector<std::vector<Cell>> &board);

public:
    explicit KakuroGenerator(int boardSize, int popSize = 50, int gens = 100,
                    double mutRate = 0.1, double crossRate = 0.7);

    std::vector<std::vector<Cell>> generateBoard(double targetFitness = 0.98,
                                                 std::chrono::seconds maxTime = std::chrono::seconds(300));

    void writeToFile(const std::string &filename, std::vector<std::vector<Cell>> &board) const;


private:
    // Tournament selection with adaptive size
    Individual tournamentSelect(const std::vector<Individual> &population, int &generation);

};


#endif //PART2_KAKUROGENERATOR_H
