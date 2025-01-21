//
// Created by Theresa on 19.01.2025.
//

#include <vector>
#include <random>
#include <algorithm>
#include <chrono>
#include "KakuroSolver.h"

class KakuroBoard {
private:
    int rows;
    int cols;

    const int minRunLength = 2;
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<> dis;
    std::uniform_int_distribution<> numDis;

public:
    std::vector<std::vector<Cell>> board;

    KakuroBoard(int r, int c) :
            rows(r),
            cols(c),
            board(r, std::vector<Cell>(c)),
            gen(rd()),
            dis(0.0, 1.0),
            numDis(1, 9) {}

    bool generateBoard() {
        const int maxAttempts = 100;
        int attempts = 0;

        while (attempts < maxAttempts) {
            // Reset board
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    board[i][j] = Cell(true); // Black cell by default
                }
            }

            placeWhiteCells();
            // check if the board is full of black fields
            bool notBlack = false;
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    if (!board[i][j].isBlack) {
                        notBlack = true;
                        break;
                    }
                }
            }
            if (!notBlack) return false;

            if (validateAndFixRuns()) {
                placeClues();
                //std::cout << " VALIDATED WITH CLUES" << std::endl;
                //printBoard();
                return true;
            }

            attempts++;
        }

        return false;
    }

    void printBoard() const {
        std::cout << "\nCurrent board state:\n";
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                const Cell &cell = board[i][j];
                std::cout << std::setw(8) << " ";  // Initial padding for all cells
                if (cell.isBlack) {
                    if (cell.downClue > 0 && cell.rightClue > 0) {
                        std::cout << std::right << cell.downClue << "\\"
                                  << std::left << std::setw(2) << cell.rightClue << "     ";
                    } else if (cell.downClue > 0) {
                        std::cout << std::right << cell.downClue << "\\" << std::left
                                  << std::setw(2) << "     ";
                    } else if (cell.rightClue > 0) {
                        std::cout << std::right << "\\" << std::left
                                  << std::setw(2) << cell.rightClue << "     ";
                    } else {
                        std::cout << "#" << "       ";
                    }
                } else {
                    std::cout << (cell.value == 0 ? "_" : std::to_string(cell.value))
                              << "       ";
                }
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }

private:


    void placeWhiteCells() {
        for (int i = 1; i < rows; ++i) {
            for (int j = 1; j < cols; ++j) {
                if (dis(gen) < 0.7) { // 70% chance of white cell
                    board[i][j] = Cell(false);
                }
            }
        }
    }

    bool validateAndFixRuns() {
        // Check horizontal runs
        for (int i = 0; i < rows; ++i) {
            int runLength = 0;
            for (int j = 0; j < cols; ++j) {
                if (!board[i][j].isBlack) {
                    runLength++;
                } else {
                    if (runLength == 1) {
                        // Fix single cells by either extending or removing
                        if (j < cols - 1 && dis(gen) < 0.5) {
                            board[i][j - 1] = Cell(true);
                        } else {
                            board[i][j - 1] = Cell(false);
                        }
                    }
                    runLength = 0;
                }
            }
            // Check end of row
            if (runLength == 1) {
                board[i][cols - 1] = Cell(true);
            }
        }

        // Check vertical runs
        for (int j = 0; j < cols; ++j) {
            int runLength = 0;
            for (int i = 0; i < rows; ++i) {
                if (!board[i][j].isBlack) {
                    runLength++;
                } else {
                    if (runLength == 1) {
                        // Fix single cells
                        if (i < rows - 1 && dis(gen) < 0.5) {
                            board[i - 1][j] = Cell(true);
                        } else {
                            board[i - 1][j] = Cell(false);
                        }
                    }
                    runLength = 0;
                }
            }
            // Check end of column
            if (runLength == 1) {
                board[rows - 1][j] = Cell(true);
            }
        }
        //std::cout << "after validate & fix runs" << std::endl;
        //printBoard();
        return validateBoard();
    }

    bool validateBoard() {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (board[i][j].isBlack) {
                    // Check horizontal run
                    if (j + 1 < cols && !board[i][j + 1].isBlack) {
                        int runLength = 0;
                        int k = j + 1;
                        while (k < cols && !board[i][k].isBlack) {
                            runLength++;
                            k++;
                        }
                        if (runLength < minRunLength) return false;
                    }

                    // Check vertical run
                    if (i + 1 < rows && !board[i + 1][j].isBlack) {
                        int runLength = 0;
                        int k = i + 1;
                        while (k < rows && !board[k][j].isBlack) {
                            runLength++;
                            k++;
                        }
                        if (runLength < minRunLength) return false;
                    }
                }
            }
        }
        //std::cout << "VALIDATED BOARD" << std::endl;
        //printBoard();
        return true;
    }

    void placeClues() {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (board[i][j].isBlack) {
                    auto [rightSum, rightCount] = calculateRightRun(i, j);
                    auto [downSum, downCount] = calculateDownRun(i, j);

                    if (rightCount >= minRunLength) board[i][j].rightClue = rightSum;
                    if (downCount >= minRunLength) board[i][j].downClue = downSum;
                }
            }
        }
    }

    std::pair<int, int> calculateRightRun(int row, int col) {
        if (col + 1 >= cols || board[row][col + 1].isBlack) return {0, 0};

        int count = 0;
        int j = col + 1;
        while (j < cols && !board[row][j].isBlack) {
            count++;
            ++j;
        }

        if (count < minRunLength) return {0, 0};

        // Generate a reasonable sum for the run length
        // For each cell, we can use 1-9, so generate a sum that makes sense
        int sum = 0;
        for (int i = 0; i < count; ++i) {
            sum += numDis(gen);
        }

        return {sum, count};
    }

    std::pair<int, int> calculateDownRun(int row, int col) {
        if (row + 1 >= rows || board[row + 1][col].isBlack) return {0, 0};

        int count = 0;
        int i = row + 1;
        while (i < rows && !board[i][col].isBlack) {
            count++;
            ++i;
        }

        if (count < minRunLength) return {0, 0};

        // Generate a reasonable sum for the run length
        int sum = 0;
        for (int i = 0; i < count; ++i) {
            sum += numDis(gen);
        }

        return {sum, count};
    }
};

class KakuroGenerator {
private:
    struct Individual {
        std::vector<std::vector<Cell>> board;
        double fitness;

        Individual(int size) : board(size, std::vector<Cell>(size, Cell(true))), fitness(0.0) {}
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


    // Check basic board constraints
    double checkBasicValidity(const std::vector<std::vector<Cell>> &board) {
        double score = 0.0;

        // Check for white cells without runs
        std::vector<std::vector<bool>> cellInRun(size, std::vector<bool>(size, false));
        int totalWhiteCells = 0;
        int whiteCellsInRuns = 0;

        // First pass: mark cells that are part of runs
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (!board[i][j].isBlack) {
                    totalWhiteCells++;
                    // Check horizontal run
                    if (j > 0 && board[i][j - 1].rightClue > 0) {
                        cellInRun[i][j] = true;
                        whiteCellsInRuns++;
                    }
                    // Check vertical run
                    if (i > 0 && board[i - 1][j].downClue > 0) {
                        cellInRun[i][j] = true;
                        whiteCellsInRuns++;
                    }
                }
            }
        }

        if (totalWhiteCells == 0) return 0.0; // All black board

        // Calculate percentage of white cells that are part of runs
        double runCoverage = static_cast<double>(whiteCellsInRuns) / totalWhiteCells;
        score += 0.3 * runCoverage;

        // Check for valid clue placement (min two cells per clue)
        int validClues = 0;
        int totalClues = 0;
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (board[i][j].isBlack) {
                    if (board[i][j].rightClue > 0) {
                        totalClues++;
                        // Count white cells to the right until we hit a black cell or boundary
                        int whiteCount = 0;
                        for (int k = j + 1; k < size && !board[i][k].isBlack; k++) {
                            whiteCount++;
                        }
                        if (whiteCount >= 2) {
                            validClues++;
                        }
                    }
                    if (board[i][j].downClue > 0) {
                        totalClues++;
                        // Count white cells downward until we hit a black cell or boundary
                        int whiteCount = 0;
                        for (int k = i + 1; k < size && !board[k][j].isBlack; k++) {
                            whiteCount++;
                        }
                        if (whiteCount >= 2) {
                            validClues++;
                        }
                    }
                }
            }
        }

        if (totalClues > 0) {
            score += 0.3 * (static_cast<double>(validClues) / totalClues);
        }

        // Check for isolated white cells
        int connectedWhiteCells = 0;
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (!board[i][j].isBlack) {
                    bool hasConnection = false;
                    // Check adjacent cells
                    if (i > 0 && !board[i - 1][j].isBlack) hasConnection = true;
                    if (i < size - 1 && !board[i + 1][j].isBlack) hasConnection = true;
                    if (j > 0 && !board[i][j - 1].isBlack) hasConnection = true;
                    if (j < size - 1 && !board[i][j + 1].isBlack) hasConnection = true;
                    if (hasConnection) connectedWhiteCells++;
                }
            }
        }

        if (totalWhiteCells > 0) {
            score += 0.3 * (static_cast<double>(connectedWhiteCells) / totalWhiteCells);
        }

        // Check for boxed-in clues
        int boxedClues = 0;
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                // Check only black cells that have clues
                if (board[i][j].isBlack && (board[i][j].downClue > 0 || board[i][j].rightClue > 0)) {
                    bool isBoxedIn = true;
                    // Check if there's at least one white cell adjacent
                    if (i > 0 && !board[i - 1][j].isBlack) isBoxedIn = false;
                    if (i < size - 1 && !board[i + 1][j].isBlack) isBoxedIn = false;
                    if (j > 0 && !board[i][j - 1].isBlack) isBoxedIn = false;
                    if (j < size - 1 && !board[i][j + 1].isBlack) isBoxedIn = false;

                    if (isBoxedIn) boxedClues++;
                }
            }
        }

        if (boxedClues > 0) {
            score += 0.1 * (static_cast<double>(boxedClues) / totalClues);
        }

        return score;
    }

    struct Run {
        int startX, startY;
        int length;
        bool horizontal;

        Run(int x, int y, int l, bool h) : startX(x), startY(y), length(l), horizontal(h) {}
    };

    double calculateFitness(Individual &individual) {
        KakuroSolver solver(size);

        // Copy board configuration to solver
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                solver.setCell(i, j,
                               individual.board[i][j].isBlack,
                               individual.board[i][j].downClue,
                               individual.board[i][j].rightClue);
            }
        }
        // First check basic validity and provide partial fitness
        double basicValidityScore = checkBasicValidity(individual.board);
        if (basicValidityScore < 0.5) {  // If board is very invalid, return partial score
            return basicValidityScore * 0.2;  // Scale down to indicate it's far from complete
        }


        // Intermediate rewards for partial improvements
        double intermediateScore = 0.0;

        // Reward for proper clue ranges (3-45 for kakuro)
        int validClueRanges = 0;
        int totalClues = 0;
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (individual.board[i][j].downClue > 0) {
                    totalClues++;
                    if (individual.board[i][j].downClue >= 3 && individual.board[i][j].downClue <= 45) {
                        validClueRanges++;
                    }
                }
                if (individual.board[i][j].rightClue > 0) {
                    totalClues++;
                    if (individual.board[i][j].rightClue >= 3 && individual.board[i][j].rightClue <= 45) {
                        validClueRanges++;
                    }
                }
            }
        }
        if (totalClues > 0) {
            intermediateScore += 0.1 * (static_cast<double>(validClueRanges) / totalClues);
        }

        // Reward for proper run lengths (2-9 cells for kakuro)
        int validRunLengths = 0;
        int totalRuns = 0;
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (individual.board[i][j].isBlack) {
                    if (individual.board[i][j].rightClue > 0) {
                        totalRuns++;
                        int runLength = 0;
                        for (int k = j + 1; k < size && !individual.board[i][k].isBlack; k++) {
                            runLength++;
                        }
                        if (runLength >= 2 && runLength <= 9) {
                            validRunLengths++;
                        }
                    }
                    if (individual.board[i][j].downClue > 0) {
                        totalRuns++;
                        int runLength = 0;
                        for (int k = i + 1; k < size && !individual.board[k][j].isBlack; k++) {
                            runLength++;
                        }
                        if (runLength >= 2 && runLength <= 9) {
                            validRunLengths++;
                        }
                    }
                }
            }
        }
        if (totalRuns > 0) {
            intermediateScore += 0.1 * (static_cast<double>(validRunLengths) / totalRuns);
        }


        // Check if board is fully valid
        if (!solver.isValidBoard()) {
            return (basicValidityScore * 0.6) + (intermediateScore * 0.4); // Return higher partial score for nearly valid boards
        }

        std::vector<std::vector<Cell>> solution;
        SolveResult result = solver.solveBoard(solution);

        double solvabilityScore = 0.0;
        double uniquenessScore = 0.0;
        double complexityScore = 0.0;
        double balanceScore = 0.0;

        // Evaluate solvability and uniqueness
        switch (result) {
            case SolveResult::UNIQUE_SOLUTION:
                solvabilityScore = 1.0;
                uniquenessScore = 1.0;
                return 1.0;
                break;
            case SolveResult::MULTIPLE_SOLUTIONS:
                solvabilityScore = 0.6;
                uniquenessScore = 0.3;
                break;
            case SolveResult::NO_SOLUTION:
                return 0.0;
        }

        // Calculate complexity score based on number of clues and their values
        totalClues = 0;
        int totalClueValue = 0;
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (individual.board[i][j].downClue > 0) {
                    totalClues++;
                    totalClueValue += individual.board[i][j].downClue;
                }
                if (individual.board[i][j].rightClue > 0) {
                    totalClues++;
                    totalClueValue += individual.board[i][j].rightClue;
                }
            }
        }

        // Prefer boards with a good number of clues (around 30% of cells)
        double optimalClueRatio = 0.3;
        double actualClueRatio = static_cast<double>(totalClues) / (size * size);
        complexityScore = 1.0 - std::abs(actualClueRatio - optimalClueRatio);

        // Calculate balance score (distribution of black vs white cells)
        int whiteCount = 0;
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (!individual.board[i][j].isBlack) {
                    whiteCount++;
                }
            }
        }
        double whiteRatio = static_cast<double>(whiteCount) / (size * size);
        double optimalWhiteRatio = 0.7; // Prefer around 70% white cells
        balanceScore = 1.0 - std::abs(whiteRatio - optimalWhiteRatio);

        // Combine scores with weights
        return WEIGHT_SOLVABILITY * solvabilityScore +
               WEIGHT_UNIQUENESS * uniquenessScore +
               WEIGHT_COMPLEXITY * complexityScore +
               WEIGHT_BALANCE * balanceScore;
    }

    Individual createRandomIndividual() {
        Individual ind(size);
        KakuroBoard newBoard(size, size);
        bool validBoard = false;

        do {
            validBoard = newBoard.generateBoard();
        } while (!validBoard);

        ind.board = newBoard.board;

        return ind;

    }

    void heavyMutation(Individual &ind) {
        std::uniform_real_distribution<> dist(0.0, 1.0);
        std::uniform_int_distribution<> clueDist(1, 45);

        // Heavy mutation - change about 30% of the board
        for (int i = 0; i < size; ++i) {  // Skip border
            for (int j = 0; j < size; ++j) {
                if (dist(rng) < 0.3) {  // 30% chance to modify each cell
                    // Flip cell type
                    ind.board[i][j].isBlack = !ind.board[i][j].isBlack;

                    if (ind.board[i][j].isBlack) {
                        // Add new clues with higher probability
                        if (j + 1 < size && !ind.board[i][j + 1].isBlack && dist(rng) < 0.8) {
                            ind.board[i][j].rightClue = clueDist(rng);
                        }
                        if (i + 1 < size && !ind.board[i + 1][j].isBlack && dist(rng) < 0.8) {
                            ind.board[i][j].downClue = clueDist(rng);
                        }
                    } else {
                        // Clear clues if cell becomes white
                        ind.board[i][j].rightClue = 0;
                        ind.board[i][j].downClue = 0;
                    }
                }
            }
        }

        // Additionally, mutate some existing clues
        for (int i = 1; i < size; ++i) {
            for (int j = 1; j < size; ++j) {
                if (ind.board[i][j].isBlack && dist(rng) < 0.4) {  // 40% chance for each black cell
                    if (ind.board[i][j].rightClue > 0) {
                        ind.board[i][j].rightClue = clueDist(rng);
                    }
                    if (ind.board[i][j].downClue > 0) {
                        ind.board[i][j].downClue = clueDist(rng);
                    }
                }
            }
        }
    }


    Individual crossover(const Individual &parent1, const Individual &parent2) {
        Individual child(size);
        std::uniform_real_distribution<> dist(0.0, 1.0);

        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                // 50% chance to inherit from each parent
                if (dist(rng) < 0.5) {
                    child.board[i][j] = parent1.board[i][j];
                } else {
                    child.board[i][j] = parent2.board[i][j];
                }
            }
        }
        return child;
    }

    // Test with different mutation - did not help
    /*
    void mutate(Individual &ind, int &generationsWithoutImprovement) {
        std::uniform_real_distribution<> dist(0.0, 1.0);
        std::uniform_int_distribution<> clueDist(3, 45);

        // Adaptive mutation rate based on fitness and stagnation
        double adaptiveMutationRate = mutationRate * (1.0 - ind.fitness);
        if (generationsWithoutImprovement > 15) {
            adaptiveMutationRate *= (1.0 + (generationsWithoutImprovement - 15) * 0.05);
        }

        // Different types of mutations with different probabilities
        for (int i = 1; i < size; ++i) {
            for (int j = 1; j < size; ++j) {
                double mutationType = dist(rng);

                if (mutationType < adaptiveMutationRate * 0.3) {  // Small clue modifications
                    if (ind.board[i][j].isBlack) {
                        // Slightly modify existing clues
                        if (ind.board[i][j].rightClue > 0) {
                            int currentClue = ind.board[i][j].rightClue;
                            // 50% chance to increment or decrement
                            if (dist(rng) < 0.5) {
                                currentClue = std::min(45, currentClue + 1);
                            } else {
                                currentClue = std::max(3, currentClue - 1);
                            }
                            ind.board[i][j].rightClue = currentClue;
                        }
                        if (ind.board[i][j].downClue > 0) {
                            int currentClue = ind.board[i][j].downClue;
                            if (dist(rng) < 0.5) {
                                currentClue = std::min(45, currentClue + 1);
                            } else {
                                currentClue = std::max(3, currentClue - 1);
                            }
                            ind.board[i][j].downClue = currentClue;
                        }
                    }
                }
                else if (mutationType < adaptiveMutationRate * 0.6) {  // Add/remove clues
                    if (ind.board[i][j].isBlack) {
                        // Toggle clues without changing cell type
                        if (j + 1 < size && !ind.board[i][j + 1].isBlack) {
                            if (ind.board[i][j].rightClue == 0 && dist(rng) < 0.3) {
                                ind.board[i][j].rightClue = clueDist(rng);
                            } else if (ind.board[i][j].rightClue > 0 && dist(rng) < 0.1) {
                                ind.board[i][j].rightClue = 0;
                            }
                        }
                        if (i + 1 < size && !ind.board[i + 1][j].isBlack) {
                            if (ind.board[i][j].downClue == 0 && dist(rng) < 0.3) {
                                ind.board[i][j].downClue = clueDist(rng);
                            } else if (ind.board[i][j].downClue > 0 && dist(rng) < 0.1) {
                                ind.board[i][j].downClue = 0;
                            }
                        }
                    }
                }
                else if (mutationType < adaptiveMutationRate) {  // Cell type changes (least frequent)
                    // Only flip cell type if it doesn't break runs
                    bool canFlip = true;
                    if (!ind.board[i][j].isBlack) {
                        // Check if this cell is part of a valid run
                        bool inHorizontalRun = (j > 0 && ind.board[i][j-1].rightClue > 0);
                        bool inVerticalRun = (i > 0 && ind.board[i-1][j].downClue > 0);
                        if (inHorizontalRun || inVerticalRun) {
                            canFlip = false;
                        }
                    }

                    if (canFlip) {
                        ind.board[i][j].isBlack = !ind.board[i][j].isBlack;
                        if (!ind.board[i][j].isBlack) {
                            ind.board[i][j].rightClue = 0;
                            ind.board[i][j].downClue = 0;
                        }
                    }
                }
            }
        }

        // Occasional repair step
        if (dist(rng) < 0.1) {
            repairInvalidRuns(ind.board);
        }
    }

    void repairInvalidRuns(std::vector<std::vector<Cell>>& board) {
        // Fix runs that are too short or too long
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (board[i][j].isBlack) {
                    // Check horizontal runs
                    if (board[i][j].rightClue > 0) {
                        int runLength = 0;
                        for (int k = j + 1; k < size && !board[i][k].isBlack; k++) {
                            runLength++;
                        }
                        if (runLength < 2 || runLength > 9) {
                            // Either remove the clue or adjust the run
                            if (runLength < 2) {
                                board[i][j].rightClue = 0;
                            } else if (runLength > 9) {
                                // Add a black cell to shorten the run
                                board[i][j + 9].isBlack = true;
                            }
                        }
                    }

                    // Check vertical runs
                    if (board[i][j].downClue > 0) {
                        int runLength = 0;
                        for (int k = i + 1; k < size && !board[k][j].isBlack; k++) {
                            runLength++;
                        }
                        if (runLength < 2 || runLength > 9) {
                            if (runLength < 2) {
                                board[i][j].downClue = 0;
                            } else if (runLength > 9) {
                                board[i + 9][j].isBlack = true;
                            }
                        }
                    }
                }
            }
        }
    }
*/

    void mutate(Individual &ind, int &generationsWithoutImprovement) {
        std::uniform_real_distribution<> dist(0.0, 1.0);
        std::uniform_int_distribution<> clueDist(1, 45);


        // Adaptive mutation rate based on fitness
        double adaptiveMutationRate = mutationRate * (1.0 - ind.fitness);

        // Adapt mutation rate, if too many generations without improvement
        if (generationsWithoutImprovement == 15) {
            //std::cout << "20 generations without best fitness improvement - adapting mutation rate old " << adaptiveMutationRate << " new ";
            adaptiveMutationRate *= (1.0 + (generationsWithoutImprovement - 20) * 0.05);
            //std::cout << " " << adaptiveMutationRate << std::endl;
        } else if (generationsWithoutImprovement == 25) {
            adaptiveMutationRate *= (1.0 + (generationsWithoutImprovement - 20) * 0.08);
        } else if (generationsWithoutImprovement == 40) {
            adaptiveMutationRate = mutationRate * (1.0 - ind.fitness);
        }

/* //Note: test to react to generationswithout improvement did not really help
        if (generationsWithoutImprovement >= 10) {
            double optimizeProb = 0.2 + (generationsWithoutImprovement * 0.01);  // Increases with stagnation
            if (dist(rng) < optimizeProb) {
                optimizeRunLengths(ind.board);
            }
        } else if (generationsWithoutImprovement >= 5) {
            constraintGuidedMutation(ind, generationsWithoutImprovement);
        } else {
*/
        // Regular mutations
        for (int i = 1; i < size; ++i) {
            for (int j = 1; j < size; ++j) {
                if (dist(rng) < adaptiveMutationRate) {
                    // Flip cell type
                    ind.board[i][j].isBlack = !ind.board[i][j].isBlack;

                    if (ind.board[i][j].isBlack) {
                        // Maybe add clues
                        if (j + 1 < size && !ind.board[i][j + 1].isBlack) {
                            ind.board[i][j].rightClue = (dist(rng) < 0.7) ? clueDist(rng) : 0;
                        }
                        if (i + 1 < size && !ind.board[i + 1][j].isBlack) {
                            ind.board[i][j].downClue = (dist(rng) < 0.7) ? clueDist(rng) : 0;
                        }
                    } else {
                        // Clear clues if cell becomes white
                        ind.board[i][j].rightClue = 0;
                        ind.board[i][j].downClue = 0;
                    }
                } else if (ind.board[i][j].isBlack) {
                    // Mutate existing clues with lower probability
                    if (dist(rng) < adaptiveMutationRate * 0.5) {
                        if (ind.board[i][j].rightClue > 0) {
                            ind.board[i][j].rightClue = clueDist(rng);
                        }
                        if (ind.board[i][j].downClue > 0) {
                            ind.board[i][j].downClue = clueDist(rng);
                        }
                    }
                }
            }
        }
    }

    void printBoard(std::vector<std::vector<Cell>> &board) const {
        std::cout << "\nCurrent board state:\n";
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                const Cell &cell = board[i][j];
                std::cout << std::setw(8) << " ";  // Initial padding for all cells
                if (cell.isBlack) {
                    if (cell.downClue > 0 && cell.rightClue > 0) {
                        std::cout << std::right << cell.downClue << "\\"
                                  << std::left << std::setw(2) << cell.rightClue << "     ";
                    } else if (cell.downClue > 0) {
                        std::cout << std::right << cell.downClue << "\\" << std::left
                                  << std::setw(2) << "     ";
                    } else if (cell.rightClue > 0) {
                        std::cout << std::right << "\\" << std::left
                                  << std::setw(2) << cell.rightClue << "     ";
                    } else {
                        std::cout << "#" << "       ";
                    }
                } else {
                    std::cout << (cell.value == 0 ? "_" : std::to_string(cell.value))
                              << "       ";
                }
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }

    //constraint guided mutation - test in order to guide towards solvable board
    bool isClueValuePossible(int clue, int length) {
        // Validate inputs
        if (length <= 0 || length > 9) {
            return false;
        }
        if (clue <= 0) {
            return false;
        }

        // Minimum possible sum for length (1+2+3+...)
        int minSum = (length * (length + 1)) / 2;
        // Maximum possible sum for length (9+8+7+...)
        int maxSum = 0;
        for (int i = 0; i < length; ++i) {
            maxSum += (9 - i);
        }

        return clue >= minSum && clue <= maxSum;
    }

    int countRunLength(const std::vector<std::vector<Cell>> &board, int startX, int startY, bool horizontal) {
        // Validate starting position
        if (startX < 0 || startX >= size || startY < 0 || startY >= size) {
            return 0;
        }

        // Check if starting cell is valid
        if (!board[startX][startY].isBlack) {
            return 0;
        }

        // Check if there's a clue
        if (horizontal && board[startX][startY].rightClue <= 0) {
            return 0;
        }
        if (!horizontal && board[startX][startY].downClue <= 0) {
            return 0;
        }

        int length = 0;
        if (horizontal) {
            int j = startY + 1;
            while (j < size && !board[startX][j].isBlack) {
                length++;
                ++j;
            }
        } else {
            int i = startX + 1;
            while (i < size && !board[i][startY].isBlack) {
                length++;
                ++i;
            }
        }

        return length;
    }

    void adjustRunOrClue(std::vector<std::vector<Cell>> &board, int x, int y, bool horizontal,
                         std::mt19937 &rng) {
        int length = countRunLength(board, x, y, horizontal);
        if (length == 0) {
            return;  // Invalid run, skip adjustment
        }

        int clue = horizontal ? board[x][y].rightClue : board[x][y].downClue;
        if (clue <= 0) {
            return;  // Invalid clue, skip adjustment
        }

        if (!isClueValuePossible(clue, length)) {
            std::uniform_real_distribution<> dist(0.0, 1.0);

            if (dist(rng) < 0.5 || length > 4) {
                // Adjust run length by adding a black cell
                if (horizontal) {
                    int splitPoint = y + std::min(2, length - 1);
                    // Ensure we're not creating a run of length 1
                    if (splitPoint < size - 1 && splitPoint > y + 1) {
                        board[x][splitPoint].isBlack = true;
                        if (dist(rng) < 0.7) {
                            // Ensure new clue is valid for remaining length
                            int remainingLength = countRunLength(board, x, splitPoint, true);
                            if (remainingLength > 0) {
                                int minSum = (remainingLength * (remainingLength + 1)) / 2;
                                int maxSum = remainingLength * 9;
                                std::uniform_int_distribution<> clueDist(minSum, maxSum);
                                board[x][splitPoint].rightClue = clueDist(rng);
                            }
                        }
                    }
                } else {
                    int splitPoint = x + std::min(2, length - 1);
                    if (splitPoint < size - 1 && splitPoint > x + 1) {
                        board[splitPoint][y].isBlack = true;
                        if (dist(rng) < 0.7) {
                            int remainingLength = countRunLength(board, splitPoint, y, false);
                            if (remainingLength > 0) {
                                int minSum = (remainingLength * (remainingLength + 1)) / 2;
                                int maxSum = remainingLength * 9;
                                std::uniform_int_distribution<> clueDist(minSum, maxSum);
                                board[splitPoint][y].downClue = clueDist(rng);
                            }
                        }
                    }
                }
            } else {
                // Adjust clue value
                int minSum = (length * (length + 1)) / 2;
                int maxSum = length * 9;  // Simplified max calculation

                if (minSum < maxSum) {  // Only adjust if valid range exists
                    std::uniform_int_distribution<> validClueDist(minSum, maxSum);
                    if (horizontal) {
                        board[x][y].rightClue = validClueDist(rng);
                    } else {
                        board[x][y].downClue = validClueDist(rng);
                    }
                }
            }
        }
    }

    void constraintGuidedMutation(Individual &ind, int generationsStuck) {
        // More aggressive fixing when stuck
        double fixProbability = std::min(0.8, 0.3 + (generationsStuck * 0.01));
        std::uniform_real_distribution<> dist(0.0, 1.0);

        // Check and fix all runs
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (ind.board[i][j].isBlack) {
                    // Check horizontal run
                    if (ind.board[i][j].rightClue > 0 && dist(rng) < fixProbability) {
                        adjustRunOrClue(ind.board, i, j, true, rng);
                    }
                    // Check vertical run
                    if (ind.board[i][j].downClue > 0 && dist(rng) < fixProbability) {
                        adjustRunOrClue(ind.board, i, j, false, rng);
                    }
                }
            }
        }
    }

    // tested in order to optimize for runs, did not help
    std::vector<Run> analyzeRunLengths(const std::vector<std::vector<Cell>> &board) {
        std::vector<Run> runs;

        // Find horizontal runs
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (board[i][j].isBlack && board[i][j].rightClue > 0) {
                    int length = 0;
                    int k = j + 1;
                    while (k < size && !board[i][k].isBlack) {
                        length++;
                        k++;
                    }
                    if (length > 0) {
                        runs.emplace_back(i, j, length, true);
                    }
                }
            }
        }

        // Find vertical runs
        for (int j = 0; j < size; ++j) {
            for (int i = 0; i < size; ++i) {
                if (board[i][j].isBlack && board[i][j].downClue > 0) {
                    int length = 0;
                    int k = i + 1;
                    while (k < size && !board[k][j].isBlack) {
                        length++;
                        k++;
                    }
                    if (length > 0) {
                        runs.emplace_back(i, j, length, false);
                    }
                }
            }
        }

        return runs;
    }

    void optimizeRunLengths(std::vector<std::vector<Cell>> &board) {
        auto runs = analyzeRunLengths(board);
        std::uniform_real_distribution<> dist(0.0, 1.0);
        std::uniform_int_distribution<> clueDist(1, 45);

        for (const Run &run: runs) {
            if (run.length > 4) {  // Too long, try to break it
                if (run.horizontal) {
                    int breakPoint = run.startY + 2 + (dist(rng) < 0.5 ? 0 : 1);
                    if (breakPoint < size - 1) {
                        // Insert black cell to break the run
                        board[run.startX][breakPoint].isBlack = true;
                        // Add new clue with 70% probability
                        if (dist(rng) < 0.7) {
                            board[run.startX][breakPoint].rightClue = clueDist(rng);
                        }
                    }
                } else {
                    int breakPoint = run.startX + 2 + (dist(rng) < 0.5 ? 0 : 1);
                    if (breakPoint < size - 1) {
                        board[breakPoint][run.startY].isBlack = true;
                        if (dist(rng) < 0.7) {
                            board[breakPoint][run.startY].downClue = clueDist(rng);
                        }
                    }
                }
            }
        }
    }

public:
    KakuroGenerator(int
                    boardSize)
            : size(boardSize),
              populationSize(80),
              mutationRate(0.1),
              targetFitness(0.83),
              rng(std::random_device{}()) {
    }

    std::vector<std::vector<Cell>> generateBoard(int maxGenerations = 2000, int maxTimeSeconds = 300) {
        auto startTime = std::chrono::steady_clock::now();

        // Initialize population
        std::vector<Individual> population;
        for (int i = 0; i < populationSize; ++i) {
            population.push_back(createRandomIndividual());
        }

        // Main evolution loop
        int generation = 0;
        Individual bestIndividual(size);
        double bestFitness = 0.0;

        while (generation < maxGenerations) {


            // Evaluate fitness
#pragma omp parallel for
            for (int i = 0; i < populationSize; ++i) {
                population[i].fitness = calculateFitness(population[i]);
            }

            // Find best individual
            auto bestIt = std::max_element(population.begin(), population.end(),
                                           [](const Individual &a, const Individual &b) {
                                               return a.fitness < b.fitness;
                                           });
            double avgFitnessOutput = 0;

            //debug output
            for (int i = 0; i < population.size(); ++i) {
                avgFitnessOutput += population[i].fitness;
            }
            avgFitnessOutput /= population.size();

            std::cout << "Generation " << generation << ": Best fitness overall = "
                      << bestFitness << " - Best fitness this run: " << bestIt->fitness << " - AVG fitness " << avgFitnessOutput << "\n";



            // found new better fitness
            if (bestIt->fitness > bestFitness) {
                bestFitness = bestIt->fitness;
                bestIndividual = *bestIt;
                std::cout << "Generation " << generation << ": Best fitness = "
                          << bestFitness << "\n";
            }

            if (bestFitness >= targetFitness) {
                std::cout << "Target fitness reached after " << generation << " generations\n";
                KakuroSolver solver(size);
                SolveResult result = solver.solveBoard(bestIndividual.board);
                if (result == SolveResult::UNIQUE_SOLUTION) break;
            }

            // Track improvement
            static int generationsWithoutImprovement = 0;
            static double lastBestFitness = 0.0;

            if (bestFitness > lastBestFitness) {
                generationsWithoutImprovement = 0;
                lastBestFitness = bestFitness;
            } else {
                generationsWithoutImprovement++;
            }


            /*
             * //Note: didn'T work as expected
             * //force population updates if very many bad boards exist
             *   int veryBadFitnessCount = 0;

            //std::cout << "current boards " << std::endl;
            for (int i = 0; i < population.size(); ++i) {
                if (population[i].fitness < 0.4) veryBadFitnessCount++;
                //debugoutput
                //printBoard(population[i].board);
            }

            if ((veryBadFitnessCount > populationSize / 2) && generationsWithoutImprovement % 5 == 0) {
                int third = (populationSize/2);
                int count = 0;
                for (int i = 0; i < populationSize; ++i) {
                    if (population[i].fitness <= 0.2 && count <= third) {
                        population[i] = createRandomIndividual();
                        count ++;
                    } else if (count >=third) {
                        break;
                    }
                }
                std::cout << "replaced " << (populationSize / 3) << " elements in population" << std::endl;
            }*/

            // Restart mechanism if stuck
            if (generationsWithoutImprovement > 40) {
                std::cout << "Stagnation detected, performing restart at generation " << generation << "\n";
                // Keep best individual
                Individual bestInd = *bestIt;

                // Reinitialize population
                for (int i = 1; i < populationSize; ++i) {
                    population[i] = createRandomIndividual();
                }
                // Keep best from previous run
                population[0] = bestInd;

                generationsWithoutImprovement = 0;
                continue;
            }


            // Create new population
            std::vector<Individual> newPopulation;

            // Calculate population diversity
            double diversityScore = 0.0;
            for (const auto& ind : population) {
                diversityScore += std::abs(ind.fitness - avgFitnessOutput);
            }
            diversityScore /= populationSize;

            // Adjust elite count based on diversity
            int baseEliteCount = populationSize / 10;
            int eliteCount = std::max(1, static_cast<int>(baseEliteCount * (1.0 - diversityScore)));

            // Sort by fitness
            std::sort(population.begin(), population.end(),
                      [](const Individual& a, const Individual& b) {
                          return a.fitness > b.fitness;
                      });

            // Elitism with diversity preservation
            for (int i = 0; i < eliteCount; ++i) {
                newPopulation.push_back(population[i]);
            }

            // Add diversity when needed
            double diversityThreshold = std::max(0.05, 0.1 * (1.0 - bestFitness));
            if (diversityScore < diversityThreshold) {
                std::cout << "Low diversity detected (score: " << diversityScore
                          << ", threshold: " << diversityThreshold << ")" << std::endl;

                int newRandomCount = populationSize / 10;  // 10% of population
                int addedCount = 0;
                std::uniform_real_distribution<> dist(0.0, 1.0);

                for (int i = 0; i < newRandomCount && newPopulation.size() < populationSize; ++i) {
                    if (dist(rng) < 0.7) {  // 70% chance for semi-random
                        // Create semi-random individual based on best
                        Individual newInd = population[0];  // Copy best
                        heavyMutation(newInd);  // Apply aggressive mutations
                        newPopulation.push_back(newInd);
                        addedCount++;
                    } else {  // 30% chance for completely random
                        newPopulation.push_back(createRandomIndividual());
                        addedCount++;
                    }
                }
                std::cout << "Added " << addedCount << " new individuals ("
                          << int(addedCount * 0.7) << " semi-random, "
                          << int(addedCount * 0.3) << " completely random)" << std::endl;
            }

            // Add some random new individuals, if moderate stagnation
            if (generationsWithoutImprovement > 25) {
                int newRandomCount = populationSize / 3; // 30 % new random individuals to freshen population
                std::cout << "25 gens without improvement: create " << newRandomCount
                          << " new individuals to freshen population" << std::endl;
                for (int i = 0; i < newRandomCount; ++i) {
                    newPopulation.push_back(createRandomIndividual());
                }
                generationsWithoutImprovement = 0;
            }

            // Fill rest with crossover and mutation
            while (newPopulation.size() < populationSize) {
                // Tournament selection
                std::uniform_int_distribution<> select(0, populationSize - 1);
                Individual &parent1 = population[select(rng)];
                Individual &parent2 = population[select(rng)];

                Individual child = crossover(parent1, parent2);
                if (generationsWithoutImprovement == 18) {
                    constraintGuidedMutation(child, generationsWithoutImprovement);
                }
                mutate(child, generationsWithoutImprovement);
                newPopulation.push_back(child);
            }

            population = std::move(newPopulation);
            generation++;
        }

        return bestIndividual.board;
    }
};