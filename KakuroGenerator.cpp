//
// Created by Theresa on 19.01.2025.
//

#include <vector>
#include <random>
#include <algorithm>
#include <chrono>
#include <map>
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
        int endX, endY;
        bool horizontal;

        Run(int x, int y, bool h) : startX(x), startY(y), horizontal(h) {}

        Run(int sx, int sy, int ex, int ey, bool h)
                : startX(sx), startY(sy), endX(ex), endY(ey), horizontal(h) {}
    };

    //helper struct to evaluate fitness
    struct FitnessComponents {
        double basicValidity;
        double runQuality;
        double clueDistribution;
        double boardBalance;
        double solvability;

        FitnessComponents() : basicValidity(0), runQuality(0),
                              clueDistribution(0), boardBalance(0),
                              solvability(0) {}
    };

    // fitness calculation helper function
    double calculateRunQualityScore(const std::vector<std::vector<Cell>> &board) {
        double score = 0.0;
        int totalRuns = 0;

        // Analyze horizontal runs
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (board[i][j].isBlack && board[i][j].rightClue > 0) {
                    totalRuns++;
                    int runLength = 0;
                    int k = j + 1;
                    while (k < size && !board[i][k].isBlack) {
                        runLength++;
                        k++;
                    }

                    // Score based on optimal run length (2-4 cells is ideal for most puzzles)
                    if (runLength >= 2 && runLength <= 4) {
                        score += 1.0;
                    } else if (runLength > 4 && runLength <= 6) {
                        score += 0.7;
                    } else if (runLength > 6) {
                        score += 0.3;
                    }

                    // Check if clue value is reasonable for the run length
                    int clue = board[i][j].rightClue;
                    int minPossible = (runLength * (runLength + 1)) / 2;  // Sum of smallest numbers
                    int maxPossible =
                            runLength * 9 - (runLength * (runLength - 1)) / 2;  // Sum of largest possible numbers

                    if (clue >= minPossible && clue <= maxPossible) {
                        score += 0.5;
                    }
                }
            }
        }

        // Analyze vertical runs (similar to horizontal)
        for (int j = 0; j < size; ++j) {
            for (int i = 0; i < size; ++i) {
                if (board[i][j].isBlack && board[i][j].downClue > 0) {
                    totalRuns++;
                    int runLength = 0;
                    int k = i + 1;
                    while (k < size && !board[k][j].isBlack) {
                        runLength++;
                        k++;
                    }

                    if (runLength >= 2 && runLength <= 4) {
                        score += 1.0;
                    } else if (runLength > 4 && runLength <= 6) {
                        score += 0.7;
                    } else if (runLength > 6) {
                        score += 0.3;
                    }

                    int clue = board[i][j].downClue;
                    int minPossible = (runLength * (runLength + 1)) / 2;
                    int maxPossible = runLength * 9 - (runLength * (runLength - 1)) / 2;

                    if (clue >= minPossible && clue <= maxPossible) {
                        score += 0.5;
                    }
                }
            }
        }

        return totalRuns > 0 ? score / (totalRuns * 1.5) : 0.0;  // Normalize score
    }

    //fitness calculation helper function
    double calculateClueDistributionScore(const std::vector<std::vector<Cell>> &board) {
        double score = 0.0;
        int totalClues = 0;
        std::map<int, int> clueFrequency;  // Track frequency of clue values

        // Collect clue statistics
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (board[i][j].isBlack) {
                    if (board[i][j].rightClue > 0) {
                        clueFrequency[board[i][j].rightClue]++;
                        totalClues++;
                    }
                    if (board[i][j].downClue > 0) {
                        clueFrequency[board[i][j].downClue]++;
                        totalClues++;
                    }
                }
            }
        }

        // Calculate distribution score
        if (totalClues > 0) {
            // Reward diversity in clue values
            double uniqueClueRatio = static_cast<double>(clueFrequency.size()) / totalClues;
            score += uniqueClueRatio * 0.5;

            // Reward reasonable clue ranges (prefer values between 3 and 45)
            int reasonableClues = 0;
            for (const auto &[clue, freq]: clueFrequency) {
                if (clue >= 3 && clue <= 45) {
                    reasonableClues += freq;
                }
            }
            score += static_cast<double>(reasonableClues) / totalClues * 0.5;
        }

        return score;
    }

    //fitness calculation helper function
    double calculateBoardBalanceScore(const std::vector<std::vector<Cell>> &board) {
        int totalCells = size * size;
        int blackCells = 0;
        int whiteCells = 0;

        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (board[i][j].isBlack) {
                    blackCells++;
                } else {
                    whiteCells++;
                }
            }
        }

        // Ideal ratio is around 30% black cells, 70% white cells
        double idealBlackRatio = 0.3;
        double actualBlackRatio = static_cast<double>(blackCells) / totalCells;

        return 1.0 - std::abs(idealBlackRatio - actualBlackRatio);
    }

    double calculateFitness(Individual &individual) {
        FitnessComponents components;

        // Calculate basic validity (existing checkBasicValidity function)
        components.basicValidity = checkBasicValidity(individual.board);

        // If board is very invalid, return early with partial score
        if (components.basicValidity < 0.3) {
            return components.basicValidity * 0.2;
        }

        // Calculate additional components
        components.runQuality = calculateRunQualityScore(individual.board);
        components.clueDistribution = calculateClueDistributionScore(individual.board);
        components.boardBalance = calculateBoardBalanceScore(individual.board);

        // Check solvability only if other components are promising
        if (components.basicValidity > 0.6 &&
            components.runQuality > 0.5 &&
            components.clueDistribution > 0.5) {

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

            std::vector<std::vector<Cell>> solution;
            SolveResult result = solver.solveBoard(solution);

            switch (result) {
                case SolveResult::UNIQUE_SOLUTION:
                    components.solvability = 1.0;
                    return 1.0; // perfect board immediate return
                case SolveResult::MULTIPLE_SOLUTIONS:
                    components.solvability = 0.3;
                    break;
                case SolveResult::NO_SOLUTION:
                    components.solvability = 0.0;
                    break;
            }
        }

        // Combine components with weights
        return components.basicValidity * 0.25 +
               components.runQuality * 0.25 +
               components.clueDistribution * 0.15 +
               components.boardBalance * 0.10 +
               components.solvability * 0.25;
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

    std::vector<Run> identifyRuns(const std::vector<std::vector<Cell>> &board) {
        std::vector<Run> runs;

        // Find horizontal runs
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (board[i][j].isBlack && board[i][j].rightClue > 0) {
                    int endJ = j + 1;
                    while (endJ < size && !board[i][endJ].isBlack) {
                        endJ++;
                    }
                    if (endJ - j > 1) {  // Run must be at least 2 cells
                        runs.emplace_back(i, j, i, endJ - 1, true);
                    }
                }
            }
        }

        // Find vertical runs
        for (int j = 0; j < size; ++j) {
            for (int i = 0; i < size; ++i) {
                if (board[i][j].isBlack && board[i][j].downClue > 0) {
                    int endI = i + 1;
                    while (endI < size && !board[endI][j].isBlack) {
                        endI++;
                    }
                    if (endI - i > 1) {  // Run must be at least 2 cells
                        runs.emplace_back(i, j, endI - 1, j, false);
                    }
                }
            }
        }

        return runs;
    }


    Individual crossover(const Individual &parent1, const Individual &parent2) {
        Individual child(size);
        std::uniform_real_distribution<> dist(0.0, 1.0);

        // First, identify valid runs in both parents
        auto runs1 = identifyRuns(parent1.board);
        auto runs2 = identifyRuns(parent2.board);

        // Initialize child with empty board (all black cells)
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                child.board[i][j] = Cell(true);  // Start with black cells
            }
        }

        // Keep track of modified regions
        std::vector<std::vector<bool>> modified(size, std::vector<bool>(size, false));

        // First phase: Copy complete runs from parents
        for (const auto &run: (dist(rng) < 0.5 ? runs1 : runs2)) {
            if (dist(rng) < 0.5) {  // 50% chance to copy a run
                const auto &sourceBoard = (dist(rng) < 0.5 ? parent1.board : parent2.board);

                // Check if we can copy this run (no conflicts with already copied runs)
                bool canCopy = true;
                if (run.horizontal) {
                    for (int j = run.startY; j <= run.endY; ++j) {
                        if (modified[run.startX][j]) {
                            canCopy = false;
                            break;
                        }
                    }
                } else {
                    for (int i = run.startX; i <= run.endX; ++i) {
                        if (modified[i][run.startY]) {
                            canCopy = false;
                            break;
                        }
                    }
                }

                // Copy the run if possible
                if (canCopy) {
                    // Copy the black cell with clues
                    if (run.horizontal) {
                        child.board[run.startX][run.startY] = sourceBoard[run.startX][run.startY];
                        modified[run.startX][run.startY] = true;

                        // Copy white cells
                        for (int j = run.startY + 1; j <= run.endY; ++j) {
                            child.board[run.startX][j] = sourceBoard[run.startX][j];
                            modified[run.startX][j] = true;
                        }
                    } else {
                        child.board[run.startX][run.startY] = sourceBoard[run.startX][run.startY];
                        modified[run.startX][run.startY] = true;

                        // Copy white cells
                        for (int i = run.startX + 1; i <= run.endX; ++i) {
                            child.board[i][run.startY] = sourceBoard[i][run.startY];
                            modified[i][run.startY] = true;
                        }
                    }
                }
            }
        }

        // add more randomness to crossover
        for (int i = 0; i < size; ++i) {  // Skip border cells
            for (int j = 0; j < size; ++j) {  // Skip border cells
                if (!modified[i][j] && dist(rng) < 0.3) {  // 30% chance of random cell
                    child.board[i][j] = Cell(dist(rng) < 0.5);  // Random cell type
                } else if (!modified[i][j]) {
                    child.board[i][j] = (dist(rng) < 0.5) ?
                                        parent1.board[i][j] : parent2.board[i][j];
                }
            }
        }

        // Second phase: Fill remaining cells from either parent
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (!modified[i][j]) {
                    // Choose from which parent to copy
                    child.board[i][j] = (dist(rng) < 0.5) ?
                                        parent1.board[i][j] : parent2.board[i][j];
                }
            }
        }

        // Third phase: Quick validation and repair
        repairInvalidCells(child.board);

        return child;
    }

    void repairInvalidCells(std::vector<std::vector<Cell>> &board) {
        // Only fix completely isolated cells
        for (int i = 1; i < size - 1; ++i) {
            for (int j = 1; j < size - 1; ++j) {
                if (!board[i][j].isBlack) {
                    bool completelyIsolated = true;
                    // Check all 8 neighboring cells
                    for (int di = -1; di <= 1; ++di) {
                        for (int dj = -1; dj <= 1; ++dj) {
                            if (di == 0 && dj == 0) continue;
                            if (!board[i + di][j + dj].isBlack) {
                                completelyIsolated = false;
                                break;
                            }
                        }
                    }
                    if (completelyIsolated) {
                        board[i][j].isBlack = true;
                    }
                }
            }
        }
    }

    // Calculate difference between two boards
    double calculateBoardDistance(const std::vector<std::vector<Cell>> &board1,
                                  const std::vector<std::vector<Cell>> &board2) {
        int differences = 0;
        int totalCells = size * size;

        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                // Compare cell types
                if (board1[i][j].isBlack != board2[i][j].isBlack) {
                    differences++;
                    continue;
                }

                // If both are black, compare clues
                if (board1[i][j].isBlack) {
                    if (board1[i][j].rightClue != board2[i][j].rightClue) differences++;
                    if (board1[i][j].downClue != board2[i][j].downClue) differences++;
                }
            }
        }

        return static_cast<double>(differences) / (totalCells * 1.5); // Normalize to [0,1]
    }

    // Calculate population diversity
    double calculatePopulationDiversity(const std::vector<Individual> &population) {
        double totalDistance = 0.0;
        int comparisons = 0;

        // Compare each individual with every other individual
        for (size_t i = 0; i < population.size(); ++i) {
            for (size_t j = i + 1; j < population.size(); ++j) {
                totalDistance += calculateBoardDistance(population[i].board, population[j].board);
                comparisons++;
            }
        }

        // Return average distance between all pairs
        return comparisons > 0 ? totalDistance / comparisons : 0.0;
    }


    // Mutation based implementation
    struct MutationStats {
        double successRate;      // Rate of successful mutations
        int totalMutations;      // Total mutation attempts
        int successfulMutations; // Mutations that improved fitness

        MutationStats() : successRate(0.5), totalMutations(0),
                          successfulMutations(0) {}
    };

    // Track mutation statistics per generation
    MutationStats mutationStats;

    // New configuration parameters
    const double MIN_MUTATION_RATE = 0.01;
    const double MAX_MUTATION_RATE = 0.3;
    const int ADAPTATION_INTERVAL = 5;  // Generations between rate adjustments

    double calculateAdaptiveMutationRate(const Individual &ind,
                                         int generationsWithoutImprovement,
                                         double populationDiversity) {
        // Start with base mutation rate modified by individual's fitness
        double rate = mutationRate * (1.0 - ind.fitness);

        // Adjust based on generations without improvement
        double stagnationFactor = 1.0 + (generationsWithoutImprovement / 20.0);
        rate *= stagnationFactor;

        // Adjust based on population diversity
        double diversityFactor = 1.0 + (1.0 - populationDiversity);
        rate *= diversityFactor;

        // Adjust based on mutation success rate
        double successFactor = 1.0;
        if (mutationStats.totalMutations > 0) {
            successFactor = 1.0 + (0.5 - mutationStats.successRate);
        }
        rate *= successFactor;

        // Ensure rate stays within bounds
        return std::min(MAX_MUTATION_RATE, std::max(MIN_MUTATION_RATE, rate));
    }

    void updateMutationStats(double oldFitness, double newFitness) {
        mutationStats.totalMutations++;
        if (newFitness > oldFitness) {
            mutationStats.successfulMutations++;
        }

        if (mutationStats.totalMutations >= 100) {  // Update rate periodically
            mutationStats.successRate = static_cast<double>(mutationStats.successfulMutations) /
                                        mutationStats.totalMutations;
            // Reset counters but keep success rate
            mutationStats.totalMutations = 0;
            mutationStats.successfulMutations = 0;
        }
    }


    void mutate(Individual &ind, int generationsWithoutImprovement, double populationDiversity) {
        double oldFitness = ind.fitness;
        std::uniform_real_distribution<> dist(0.0, 1.0);
        std::uniform_int_distribution<> clueDist(1, 45);

        // Calculate adaptive rate
        double adaptiveMutationRate = calculateAdaptiveMutationRate(
                ind, generationsWithoutImprovement, populationDiversity);

        // Apply different mutation strategies based on board state
        bool boardModified = false;

        // Strategy 1: Fix isolated cells
        if (dist(rng) < adaptiveMutationRate * 1.5) {
            boardModified |= fixIsolatedCells(ind);
        }

        // Strategy 2: Optimize run lengths
        if (dist(rng) < adaptiveMutationRate * 1.2) {
            boardModified |= optimizeRunLengths(ind);
        }

        // Strategy 3: Regular cell mutations with adaptive rate
        for (int i = 1; i < size; ++i) {
            for (int j = 1; j < size; ++j) {
                if (dist(rng) < adaptiveMutationRate) {
                    boardModified = true;
                    // Flip cell type
                    ind.board[i][j].isBlack = !ind.board[i][j].isBlack;

                    if (ind.board[i][j].isBlack) {
                        // Add clues with higher probability for black cells
                        if (j + 1 < size && !ind.board[i][j + 1].isBlack) {
                            ind.board[i][j].rightClue = (dist(rng) < 0.7) ? clueDist(rng) : 0;
                        }
                        if (i + 1 < size && !ind.board[i + 1][j].isBlack) {
                            ind.board[i][j].downClue = (dist(rng) < 0.7) ? clueDist(rng) : 0;
                        }
                    } else {
                        // Clear clues for white cells
                        ind.board[i][j].rightClue = 0;
                        ind.board[i][j].downClue = 0;
                    }
                }
            }
        }

        // Strategy 4: Mutate existing clues
        if (boardModified) {
            for (int i = 1; i < size; ++i) {
                for (int j = 1; j < size; ++j) {
                    if (ind.board[i][j].isBlack && dist(rng) < adaptiveMutationRate * 0.5) {
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

        // Update mutation statistics
        if (boardModified) {
            double newFitness = calculateFitness(ind);
            updateMutationStats(oldFitness, newFitness);
        }
    }

    bool fixIsolatedCells(Individual &ind) {
        bool modified = false;
        for (int i = 1; i < size - 1; ++i) {
            for (int j = 1; j < size - 1; ++j) {
                if (!ind.board[i][j].isBlack) {
                    // Check if cell is isolated
                    bool hasConnection = false;
                    if (!ind.board[i - 1][j].isBlack || !ind.board[i + 1][j].isBlack ||
                        !ind.board[i][j - 1].isBlack || !ind.board[i][j + 1].isBlack) {
                        hasConnection = true;
                    }

                    if (!hasConnection) {
                        // Convert to black cell
                        ind.board[i][j].isBlack = true;
                        modified = true;
                    }
                }
            }
        }
        return modified;
    }

    bool optimizeRunLengths(Individual &ind) {
        bool modified = false;
        std::uniform_real_distribution<> dist(0.0, 1.0);
        std::uniform_int_distribution<> clueDist(1, 45);

        // Check horizontal runs
        for (int i = 0; i < size; ++i) {
            int runLength = 0;
            int startJ = -1;

            for (int j = 0; j < size; ++j) {
                if (!ind.board[i][j].isBlack) {
                    if (runLength == 0) startJ = j;
                    runLength++;
                } else {
                    if (runLength > 6) {  // Split long runs
                        int splitPoint = startJ + runLength / 2;
                        ind.board[i][splitPoint].isBlack = true;
                        if (dist(rng) < 0.7) {
                            ind.board[i][splitPoint].rightClue = clueDist(rng);
                        }
                        modified = true;
                    }
                    runLength = 0;
                }
            }
        }

        // Check vertical runs (similar to horizontal)
        for (int j = 0; j < size; ++j) {
            int runLength = 0;
            int startI = -1;

            for (int i = 0; i < size; ++i) {
                if (!ind.board[i][j].isBlack) {
                    if (runLength == 0) startI = i;
                    runLength++;
                } else {
                    if (runLength > 6) {
                        int splitPoint = startI + runLength / 2;
                        ind.board[splitPoint][j].isBlack = true;
                        if (dist(rng) < 0.7) {
                            ind.board[splitPoint][j].downClue = clueDist(rng);
                        }
                        modified = true;
                    }
                    runLength = 0;
                }
            }
        }

        return modified;
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

    // Island model for evolving
    struct Island {
        std::vector<Individual> population;
        double bestFitness;
        Individual bestIndividual;
        double mutationRate;  // Island-specific mutation rate
        double crossoverRate; // Island-specific crossover rate
        int generationsWithoutImprovement;

        Island(int size, int popSize) :
                population(popSize, Individual(size)),
                bestFitness(0.0),
                bestIndividual(size),
                mutationRate(0.1),
                crossoverRate(0.7),
                generationsWithoutImprovement(0) {}
    };

    // Configuration for island model
    static const int NUM_ISLANDS = 6;
    static const int MIGRATION_INTERVAL = 25;
    static const int MIGRATION_SIZE = 2;
    static const int MIN_ISLAND_SIZE = 15;
    static const int OPTIMAL_ISLAND_SIZE = 20;
    const double MIGRATION_SELECTION_PRESSURE = 0.8;

    std::vector<Island> islands;

    double calculateAverageIslandFitness(const Island &island) {
        if (island.population.empty()) return 0.0;
        double sum = 0.0;
        for (const auto &ind: island.population) {
            sum += ind.fitness;
        }
        return sum / island.population.size();
    }


    void initializeWithDensity(Island &island, double blackCellProbability) {
        for (auto &individual: island.population) {
            std::uniform_real_distribution<> dist(0.0, 1.0);

            // Initialize board with specified density
            for (int i = 1; i < size - 1; ++i) {
                for (int j = 1; j < size - 1; ++j) {
                    if (dist(rng) < blackCellProbability) {
                        individual.board[i][j] = Cell(true);
                    } else {
                        individual.board[i][j] = Cell(false);
                    }
                }
            }

            // Calculate initial fitness
            individual.fitness = calculateFitness(individual);

            // Update island's best if needed
            if (individual.fitness > island.bestFitness) {
                island.bestFitness = individual.fitness;
                island.bestIndividual = individual;
            }
        }
    }

    void initializeIslands() {
        islands.clear();
        const int islandSize = populationSize / NUM_ISLANDS;

        for (int i = 0; i < NUM_ISLANDS; ++i) {
            islands.emplace_back(size, islandSize);

            // Set island-specific parameters
            islands[i].mutationRate = mutationRate * (0.5 + (i / (NUM_ISLANDS - 1.0)));
            islands[i].crossoverRate = 0.6 + (0.3 * i / (NUM_ISLANDS - 1.0));

            // Different initialization strategy per island
            switch (i % 3) {
                case 0:  // Dense black cells
                    initializeWithDensity(islands[i], 0.4);
                    break;
                case 1:  // Sparse black cells
                    initializeWithDensity(islands[i], 0.2);
                    break;
                case 2:  // Random density
                    initializeWithDensity(islands[i], 0.3);
                    break;
            }
        }
    }

    void balanceIslandSize(Island &island) {
        while (island.population.size() < MIN_ISLAND_SIZE) {
            Individual newInd = island.bestIndividual;
            heavyMutation(newInd);
            newInd.fitness = calculateFitness(newInd);
            island.population.push_back(newInd);
        }

        while (island.population.size() > OPTIMAL_ISLAND_SIZE) {
            auto worst = std::min_element(island.population.begin(),
                                          island.population.end(),
                                          [](const Individual &a, const Individual &b) {
                                              return a.fitness < b.fitness;
                                          });
            island.population.erase(worst);
        }
    }

    void performMigration() {
        // Check if migration is worthwhile
        std::vector<double> avgFitnesses(NUM_ISLANDS);
        for (int i = 0; i < NUM_ISLANDS; ++i) {
            avgFitnesses[i] = calculateAverageIslandFitness(islands[i]);
        }

        double maxDiff = 0.0;
        for (int i = 0; i < NUM_ISLANDS; ++i) {
            for (int j = i + 1; j < NUM_ISLANDS; ++j) {
                maxDiff = std::max(maxDiff, std::abs(avgFitnesses[i] - avgFitnesses[j]));
            }
        }

        if (maxDiff < 0.1) {
            return;  // Skip migration if islands are too similar
        }

        std::uniform_real_distribution<> dist(0.0, 1.0);

        // Create copies of migrants before migration
        std::vector<std::vector<Individual>> migrants(NUM_ISLANDS);

        // First select migrants from each island
        for (int sourceIsland = 0; sourceIsland < NUM_ISLANDS; ++sourceIsland) {
            auto &source = islands[sourceIsland].population;

            // Sort by fitness
            std::sort(source.begin(), source.end(),
                      [](const Individual &a, const Individual &b) {
                          return a.fitness > b.fitness;
                      });

            // Select individuals to migrate
            for (int m = 0; m < MIGRATION_SIZE; ++m) {
                int sourceIndex;
                if (dist(rng) < MIGRATION_SELECTION_PRESSURE) {
                    sourceIndex = m;  // Select from top individuals
                } else {
                    sourceIndex = std::uniform_int_distribution<>(
                            0, source.size() - 1)(rng);
                }
                migrants[sourceIsland].push_back(source[sourceIndex]);
            }
        }

        // Then perform migration using copies
        for (int sourceIsland = 0; sourceIsland < NUM_ISLANDS; ++sourceIsland) {
            int destIsland = (sourceIsland + 1) % NUM_ISLANDS;  // Ring topology
            auto &dest = islands[destIsland].population;

            // Replace worst individuals in destination
            std::sort(dest.begin(), dest.end(),
                      [](const Individual &a, const Individual &b) {
                          return a.fitness < b.fitness;  // Sort ascending for replacement
                      });

            for (int m = 0; m < MIGRATION_SIZE; ++m) {
                dest[m] = migrants[sourceIsland][m];
            }
        }
    }


    void evolveIsland(Island &island, int generation) {
        // Calculate population diversity for this island
        double populationDiversity = calculatePopulationDiversity(island.population);

        // Evaluate fitness
#pragma omp parallel for
        for (int i = 0; i < island.population.size(); ++i) {
            double newFitness = calculateFitness(island.population[i]);
            island.population[i].fitness = newFitness;

#pragma omp critical
            {
                if (newFitness > island.bestFitness) {
                    island.bestFitness = newFitness;
                    island.bestIndividual = island.population[i];
                    island.generationsWithoutImprovement = 0;
                }
            }
        }

        // Create new population
        std::vector<Individual> newPopulation;

        // Elitism
        int eliteCount = std::max(1, static_cast<int>(
                island.population.size() * 0.1 * (1.0 - populationDiversity)));

        std::sort(island.population.begin(), island.population.end(),
                  [](const Individual &a, const Individual &b) {
                      return a.fitness > b.fitness;
                  });

        for (int i = 0; i < eliteCount; ++i) {
            newPopulation.push_back(island.population[i]);
        }

        // Add diversity when needed - TODO: maybe adapt to islands, needs to be tested
        /* double diversityThreshold = std::max(0.4, 0.1 * (1.0 - island.bestFitness));
         if (populationDiversity < diversityThreshold) {
             std::cout << "Low diversity detected (score: " << populationDiversity
                       << ", threshold: " << diversityThreshold << ")" << std::endl;

             int newRandomCount = populationSize / 10;  //
             int addedCount = 0;
             std::uniform_real_distribution<> dist(0.0, 1.0);

             for (int i = 0; i < newRandomCount && newPopulation.size() < populationSize; ++i) {
                 if (dist(rng) < 0.7) {  // 70% chance for semi-random
                     // Create semi-random individual based on best
                     Individual newInd = island.bestIndividual;  // Copy best
                     heavyMutation(newInd);  // Apply aggressive mutations
                     newPopulation.push_back(newInd);
                     addedCount++;
                 } else {  // 30% chance for completely random
                     newPopulation.push_back(createRandomIndividual());
                     addedCount++;
                 }
             }
         }*/

        // Fill rest with crossover and mutation
        std::uniform_real_distribution<> dist(0.0, 1.0);
        while (newPopulation.size() < island.population.size()) {
            // Tournament selection
            std::uniform_int_distribution<> select(0, island.population.size() - 1);
            Individual &parent1 = island.population[select(rng)];
            Individual &parent2 = island.population[select(rng)];

            Individual child(size);
            if (dist(rng) < island.crossoverRate) {
                child = crossover(parent1, parent2);
            } else {
                child = dist(rng) < 0.5 ? parent1 : parent2;
            }

            if (dist(rng) < island.mutationRate) {
                mutate(child, island.generationsWithoutImprovement, populationDiversity);
            }

            newPopulation.push_back(child);
        }

        island.population = std::move(newPopulation);
        island.generationsWithoutImprovement++;

        // Balance population size
        balanceIslandSize(island);
    }


    void reinitializeIslandIfStuck(Island &island, int maxStagnantGenerations = 30) {
        if (island.generationsWithoutImprovement > maxStagnantGenerations) {
            std::cout << "Reinitializing stagnant island (no improvement for "
                      << maxStagnantGenerations << " generations)\n";

            // Keep best individual
            Individual best = island.bestIndividual;

            // Reinitialize with random individuals
            for (auto &ind: island.population) {
                ind = createRandomIndividual();
                ind.fitness = calculateFitness(ind);
            }

            // Add back best individual
            island.population[0] = best;
            island.generationsWithoutImprovement = 0;
        }
    }


public:
    KakuroGenerator(int
                    boardSize)
            : size(boardSize),
              populationSize(120),
              mutationRate(0.1),
              targetFitness(0.83),
              rng(std::random_device{}()) {
    }


    std::vector<std::vector<Cell>> generateBoard(int maxGenerations = 2000, int maxTimeSeconds = 300) {
        auto startTime = std::chrono::steady_clock::now();

        // Initialize islands
        initializeIslands();

        // Main evolution loop
        int generation = 0;
        Individual globalBest(size);
        double globalBestFitness = 0.0;

        while (generation < maxGenerations) {
            // Evolve each island
#pragma omp parallel for
            for (int i = 0; i < NUM_ISLANDS; ++i) {
                if (islands[i].population.empty()) {
                    std::cout << "Warning: Island " << i << " has empty population!\n";
                    continue;
                }
                evolveIsland(islands[i], generation);
                reinitializeIslandIfStuck(islands[i]);
            }

            // Perform migration if needed
            if (generation % MIGRATION_INTERVAL == 0 && generation > 0) {
                performMigration();
            }


            //Migration debug output
            for (const auto &island: islands) {
                double islandPopSize = island.population.size();
                double islandAvgFitness = 0.0;
                for (const auto &ind: island.population) {
                    islandAvgFitness += ind.fitness;
                }

                islandAvgFitness /= islandPopSize;
                //std::cout << "Island avg fitness: " << islandAvgFitness<<"\n";
                //<< " - Island pop size: " << std::to_string(islandPopSize) << "\n";
            }

            // Update global best
            for (const auto &island: islands) {
                if (island.bestFitness > globalBestFitness) {
                    globalBestFitness = island.bestFitness;
                    globalBest = island.bestIndividual;
                    std::cout << "Generation " << generation << ": New best fitness = "
                              << globalBestFitness << " (Island "
                              << &island - &islands[0] << ")\n";
                }
            }

            // Calculate and output statistics
            double avgFitness = 0.0;
            int totalPopulation = 0;
            for (const auto &island: islands) {
                for (const auto &ind: island.population) {
                    avgFitness += ind.fitness;
                    totalPopulation++;
                }
            }
            avgFitness /= totalPopulation;

            // Calculate global diversity
            std::vector<Individual> combinedPop;
            for (const auto &island: islands) {
                combinedPop.insert(combinedPop.end(),
                                   island.population.begin(),
                                   island.population.end());
            }
            double globalDiversity = calculatePopulationDiversity(combinedPop);

            std::cout << "Generation " << generation
                      << ": Best fitness overall = " << globalBestFitness
                      << " - AVG fitness " << avgFitness
                      << " - Diversity " << globalDiversity << "\n";

            if (globalBestFitness >= targetFitness) {
                std::cout << "Target fitness reached after " << generation << " generations\n";
                printBoard(globalBest.board);
                KakuroSolver solver(size);
                SolveResult result = solver.solveBoard(globalBest.board);
                if (result == SolveResult::UNIQUE_SOLUTION) break;
                globalBestFitness = 0;
                globalBest = Individual(size);
            }

            generation++;
        }

        return globalBest.board;
    }
};