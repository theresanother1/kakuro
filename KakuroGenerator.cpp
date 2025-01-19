//
// Created by Theresa on 19.01.2025.
//

#include "KakuroGenerator.h"

// Helper function to generate a random board layout
std::vector<std::vector<Cell>> KakuroGenerator::generateRandomBoard() {

    std::vector<std::vector<Cell>> board(size, std::vector<Cell>(size));

    // Always make the top-left cell black
    board[0][0] = Cell(true);

    // Generate random black/white cells
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i == 0 || j == 0) {
                board[i][j] = Cell(true);  // First row and column are always black
            } else {
                // Random chance of being black (30%)
                board[i][j] = Cell(std::uniform_real_distribution<>(0, 1)(rng) < 0.3);
            }
        }
    }

    // Generate random clues for black cells
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (board[i][j].isBlack) {
                // Check if we can place right clue
                int rightLength = 0;
                for (int k = j + 1; k < size && !board[i][k].isBlack; k++) {
                    rightLength++;
                }
                if (rightLength > 0) {
                    board[i][j].rightClue = generateValidClue(rightLength);
                }

                // Check if we can place down clue
                int downLength = 0;
                for (int k = i + 1; k < size && !board[k][j].isBlack; k++) {
                    downLength++;
                }
                if (downLength > 0) {
                    board[i][j].downClue = generateValidClue(downLength);
                }
            }
        }
    }

    return board;
}

// Generate a valid clue for a given run length
int KakuroGenerator::generateValidClue(int length) {
    // Calculate minimum and maximum possible sums for this length
    int minSum = (length * (length + 1)) / 2;  // Sum of smallest possible numbers
    int maxSum = (length * (19 - length)) / 2;  // Sum of largest possible numbers

    return std::uniform_int_distribution<>(minSum, maxSum)(rng);
}


// Calculate fitness of a board
double KakuroGenerator::calculateFitness(const std::vector<std::vector<Cell>> &board, int &generation) {
    double fitness = 0.0;
    std::cout << "Current generation: " << generation << std::endl;
    // Create a temporary solver to check the board
    KakuroSolver solver(size);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            solver.setCell(i, j, board[i][j].isBlack,
                           board[i][j].downClue, board[i][j].rightClue);
        }
    }

    // Criteria for fitness:
    // 1. Connected white cells
    fitness += evaluateConnectivity(board) * 0.3;

    // 2. Valid clue distribution
    fitness += evaluateClueDistribution(board) * 0.2;

    // 3. Solvability
    std::vector<std::vector<Cell>> solution;
    bool solvableUnique = solver.solveBoard(solution);
    if (solvableUnique) {
        fitness += 1;
    }

    return fitness;
}

// Evaluate connectivity of white cells
double KakuroGenerator::evaluateConnectivity(const std::vector<std::vector<Cell>> &board) {
    std::vector<std::vector<bool>> visited(size, std::vector<bool>(size, false));
    int components = 0;
    int totalWhiteCells = 0;

    // Count connected components of white cells
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (!board[i][j].isBlack) {
                totalWhiteCells++;
                if (!visited[i][j]) {
                    components++;
                    dfs(board, visited, i, j);
                }
            }
        }
    }

    // Perfect score if all white cells are connected
    return components == 1 ? 1.0 : 1.0 / components;
}

// DFS helper for connectivity check
void KakuroGenerator::dfs(const std::vector<std::vector<Cell>> &board,
         std::vector<std::vector<bool>> &visited, int i, int j) {
    if (i < 0 || i >= size || j < 0 || j >= size ||
        board[i][j].isBlack || visited[i][j]) {
        return;
    }

    visited[i][j] = true;
    dfs(board, visited, i + 1, j);
    dfs(board, visited, i - 1, j);
    dfs(board, visited, i, j + 1);
    dfs(board, visited, i, j - 1);
}

// Evaluate distribution of clues
double KakuroGenerator::evaluateClueDistribution(const std::vector<std::vector<Cell>> &board) {
    int totalClues = 0;
    int validClues = 0;

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (board[i][j].isBlack) {
                if (board[i][j].rightClue > 0) {
                    totalClues++;
                    // Check if right clue is valid
                    int length = 0;
                    for (int k = j + 1; k < size && !board[i][k].isBlack; k++) {
                        length++;
                    }
                    if (length > 0 && isValidClue(board[i][j].rightClue, length)) {
                        validClues++;
                    }
                }
                if (board[i][j].downClue > 0) {
                    totalClues++;
                    // Check if down clue is valid
                    int length = 0;
                    for (int k = i + 1; k < size && !board[k][j].isBlack; k++) {
                        length++;
                    }
                    if (length > 0 && isValidClue(board[i][j].downClue, length)) {
                        validClues++;
                    }
                }
            }
        }
    }

    return totalClues > 0 ? static_cast<double>(validClues) / totalClues : 0.0;
}

// Check if a clue is valid for a given length
bool KakuroGenerator::isValidClue(int sum, int length) {
    int minSum = (length * (length + 1)) / 2;
    int maxSum = (length * (19 - length)) / 2;
    return sum >= minSum && sum <= maxSum;
}

// Crossover two boards to create a new one
std::vector<std::vector<Cell>> KakuroGenerator::crossover(const std::vector<std::vector<Cell>> &parent1,
                                                          const std::vector<std::vector<Cell>> &parent2) {
    std::vector<std::vector<Cell>> child(size, std::vector<Cell>(size));

    // Random crossover points for both rows and columns
    int rowCrossPoint = std::uniform_int_distribution<>(1, size - 1)(rng);
    int colCrossPoint = std::uniform_int_distribution<>(1, size - 1)(rng);

    // Copy layout from parents using both row and column crossover points
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i < rowCrossPoint && j < colCrossPoint) {
                child[i][j] = parent1[i][j];  // Top-left quadrant from parent1
            } else if (i < rowCrossPoint && j >= colCrossPoint) {
                child[i][j] = parent2[i][j];  // Top-right quadrant from parent2
            } else if (i >= rowCrossPoint && j < colCrossPoint) {
                child[i][j] = parent2[i][j];  // Bottom-left quadrant from parent2
            } else {
                child[i][j] = parent1[i][j];  // Bottom-right quadrant from parent1
            }
        }
    }

    return child;
}

// Mutate a board
void KakuroGenerator::mutate(std::vector<std::vector<Cell>> &board) {
    for (int i = 1; i < size; i++) {  // Skip first row
        for (int j = 1; j < size; j++) {  // Skip first column
            if (std::uniform_real_distribution<>(0, 1)(rng) < mutationRate) {
                // Flip cell color
                board[i][j].isBlack = !board[i][j].isBlack;

                if (board[i][j].isBlack) {
                    // Generate new clues if needed
                    int rightLength = 0;
                    for (int k = j + 1; k < size && !board[i][k].isBlack; k++) {
                        rightLength++;
                    }
                    if (rightLength > 0) {
                        board[i][j].rightClue = generateValidClue(rightLength);
                    }

                    int downLength = 0;
                    for (int k = i + 1; k < size && !board[k][j].isBlack; k++) {
                        downLength++;
                    }
                    if (downLength > 0) {
                        board[i][j].downClue = generateValidClue(downLength);
                    }
                } else {
                    // Clear clues
                    board[i][j].rightClue = 0;
                    board[i][j].downClue = 0;
                }
            }
        }
    }
}


KakuroGenerator::KakuroGenerator(int boardSize, int popSize, int gens,
                double mutRate, double crossRate)
        : size(boardSize), populationSize(popSize), generations(gens),
          mutationRate(mutRate), crossoverRate(crossRate),
          rng(std::random_device{}()) {}

std::vector<std::vector<Cell>> KakuroGenerator::generateBoard(double targetFitness,
                                             std::chrono::seconds maxTime) {
    // Initialize population
    std::vector<Individual> population(populationSize, Individual(size));
    for (auto &ind: population) {
        ind.board = generateRandomBoard();
        int generation_0 = 0;
        ind.fitness = calculateFitness(ind.board, generation_0);
    }

    // Track generation count and start time
    int generation = 0;
    auto startTime = std::chrono::steady_clock::now();
    double bestFitnessAchieved = 0.0;

    // Main evolutionary loop - run until we find a good solution or timeout
    while (true) {
        std::vector<Individual> newPopulation;

        // Elitism: keep best individual
        auto bestIt = std::max_element(population.begin(), population.end(),
                                       [](const Individual &a, const Individual &b) {
                                           return a.fitness < b.fitness;
                                       });
        newPopulation.push_back(*bestIt);

        // Check if we've found a good enough solution
        bestFitnessAchieved = bestIt->fitness;
        if (bestFitnessAchieved >= targetFitness) {
            std::cout << "Target fitness achieved after " << generation
                      << " generations with fitness " << bestFitnessAchieved << std::endl;
            return bestIt->board;
        }

        /*// Check timeout
        auto currentTime = std::chrono::steady_clock::now();
        auto elapsedTime = std::chrono::duration_cast<std::chrono::seconds>(
                currentTime - startTime);
        if (elapsedTime > maxTime) {
            std::cout << "Time limit reached after " << generation
                      << " generations. Best fitness achieved: "
                      << bestFitnessAchieved << std::endl;
            return bestIt->board;
        }*/

        // Generate new individuals
        while (newPopulation.size() < populationSize) {
            // Tournament selection
            auto parent1 = tournamentSelect(population, generation);
            auto parent2 = tournamentSelect(population, generation);

            Individual child(size);
            if (std::uniform_real_distribution<>(0, 1)(rng) < crossoverRate) {
                child.board = crossover(parent1.board, parent2.board);
            } else {
                child.board = parent1.board;
            }

            mutate(child.board);
            child.fitness = calculateFitness(child.board, generation);

            // fast return, if uniquely solvable board is found
            if (child.fitness > targetFitness) return child.board;
            newPopulation.push_back(child);
        }
        std::cout << "created new population in generation " << generation << std::endl;


        // Replace 5% of population with new random boards (excluding the best individual)
        int numToReplace = std::max(1, static_cast<int>(populationSize * 0.05));
        for (int i = 0; i < numToReplace; i++) {
            int idx = std::uniform_int_distribution<>(1, populationSize - 1)(rng); // Skip index 0 (best individual)
            Individual newInd(size);
            newInd.board = generateRandomBoard();
            newInd.fitness = calculateFitness(newInd.board, generation);
            population[idx] = newInd;
        }

        population = std::move(newPopulation);

        // Print progress with detailed fitness breakdown
        if (generation % 10 == 0) {
            // Get detailed fitness components for best individual
            double connectivityScore = evaluateConnectivity(bestIt->board);
            double clueScore = evaluateClueDistribution(bestIt->board);
            std::vector<std::vector<Cell>> solution;
            bool isSolvable = KakuroSolver(size).solveBoard(solution);

            std::cout << "\nDetailed fitness breakdown for best individual:"
                      << "\n - Connectivity (30%): " << (connectivityScore * 0.3)
                      << "\n - Clue Distribution (20%): " << (clueScore * 0.2)
                      << "\n - Solvability (50%): " << (isSolvable ? 0.5 : 0.0)
                      << "\n - Total: " << bestFitnessAchieved << std::endl;
            std::cout << "Generation " << generation << ": Best fitness = "
                      << bestFitnessAchieved << std::endl;

            // Additional diagnostic information
            auto currentTime = std::chrono::steady_clock::now();
            auto elapsedTime = std::chrono::duration_cast<std::chrono::seconds>(
                    currentTime - startTime);
            std::cout << "Time elapsed: " << elapsedTime.count() << "s" << std::endl;

            // Print population diversity statistics
            double avgFitness = 0.0;
            double minFitness = 0.98;
            for (const auto &ind: population) {
                avgFitness += ind.fitness;
                minFitness = std::min(minFitness, ind.fitness);
            }
            avgFitness /= population.size();
            std::cout << "Population stats - Avg fitness: " << avgFitness
                      << ", Min fitness: " << minFitness << std::endl;
        }

        generation++;

        // Adjust mutation rate based on population diversity
        if (generation % 50 == 0) {
            double avgFitness = 0.0;
            for (const auto &ind: population) {
                avgFitness += ind.fitness;
            }
            avgFitness /= population.size();

            // If population is converging, increase mutation rate
            if (bestFitnessAchieved - avgFitness < 0.1) {
                mutationRate = std::min(mutationRate * 1.5, 0.5);
                std::cout << "Increasing mutation rate to " << mutationRate << std::endl;
            } else {
                mutationRate = std::max(mutationRate * 0.9, 0.05);
            }
        }
    }
}


// Tournament selection with adaptive size
KakuroGenerator::Individual KakuroGenerator::tournamentSelect(const std::vector<Individual> &population, int &generation) {
    // Increase tournament size as generations progress to increase selection pressure
    const int baseTournamentSize = 3;
    const int maxTournamentSize = 7;
    int tournamentSize = std::min(baseTournamentSize + (generation / 100), maxTournamentSize);
    Individual best = population[std::uniform_int_distribution<>(
            0, population.size() - 1)(rng)];

    for (int i = 1; i < tournamentSize; i++) {
        Individual contestant = population[std::uniform_int_distribution<>(
                0, population.size() - 1)(rng)];
        if (contestant.fitness > best.fitness) {
            best = contestant;
        }
    }

    return best;
}

void KakuroGenerator::writeToFile(const std::string &filename, std::vector<std::vector<Cell>> &board) const {
    std::ofstream outFile(filename);
    if (!outFile) {
        throw std::runtime_error("Could not open file for writing: " + filename);
    }

    // Write dimensions
    outFile << size << " " << size << std::endl;

    // Write board
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            const Cell &cell = board[i][j];
            outFile << std::setw(8) << " ";  // Initial padding
            if (cell.isBlack) {
                if (cell.downClue > 0 && cell.rightClue > 0) {
                    outFile << std::right << cell.downClue << "\\"
                            << std::left << std::setw(2) << cell.rightClue << "     ";
                } else if (cell.downClue > 0) {
                    outFile << std::right << cell.downClue << "\\" << std::left
                            << std::setw(2)  << "     ";
                } else if (cell.rightClue > 0) {
                    outFile << std::right << "\\" << std::left
                            << std::setw(2) << cell.rightClue << "     ";
                } else {
                    outFile << "#" << "       ";
                }
            } else {
                outFile << (cell.value == 0 ? "_" : std::to_string(cell.value))
                        << "       ";
            }
        }
        outFile << std::endl;
    }
    outFile.close();
    std::cout << "Board successfully written to " << filename << std::endl;
}


