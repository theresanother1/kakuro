#include "KakuroSolver.h"
#include <queue>
#include <atomic>

void KakuroSolver::precomputeSumCombinations() {
    // Initialize min/max sums for quick lookup
    std::vector<int> minSum(10, 0);
    std::vector<int> maxSum(10, 0);

    for (int len = 1; len <= 9; len++) {
        // Minimum sum: 1+2+3+...+len
        minSum[len] = (len * (len + 1)) / 2;
        // Maximum sum: 9+8+7+...+(9-len+1)
        maxSum[len] = (len * (19 - len)) / 2;
    }

    // Generate all valid combinations for each sum and length
    for (int length = 1; length <= 9; ++length) {
        for (int sum = minSum[length]; sum <= maxSum[length]; ++sum) {
            generateOptimizedCombinations(sum, length, sumCombinations[sum * 100 + length]);
        }
    }
}


void KakuroSolver::generateOptimizedCombinations(int targetSum, int length, std::vector<std::vector<int>> &result) {
    std::vector<int> current(length);
    std::vector<bool> used(10, false);  // Numbers 1-9

    // Helper function to generate permutations

    std::function<void(int, int)> generatePermutations = [&](int pos, int remainingSum) {
        if (pos == length) {
            if (remainingSum == 0) {
                result.push_back(current);
            }
            return;
        }

        for (int num = 1; num <= 9; num++) {
            if (!used[num] && num <= remainingSum) {
                current[pos] = num;
                used[num] = true;
                generatePermutations(pos + 1, remainingSum - num);
                used[num] = false;
            }
        }
    };

    generatePermutations(0, targetSum);
}

void KakuroSolver::printBoard() const {
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

bool KakuroSolver::isValidInRun(const Run &run, const std::vector<int> &values) {
    std::set<int> used;
    int currentSum = 0;
    size_t valueIndex = 0;

    for (size_t i = 0; i < run.cells.size() && valueIndex < values.size(); i++) {
        const auto &[x, y] = run.cells[i];
        if (!board[x][y].isBlack) {
            int value = values[valueIndex++];
            if (used.count(value)) return false;
            used.insert(value);
            currentSum += value;
        }
    }

    return currentSum == run.sum;
}

std::vector<std::vector<int>> KakuroSolver::getPossibleValues(const Run& run) {
    int key = run.sum * 100 + run.length;
    std::vector<std::vector<int>> validCombinations;

#pragma omp parallel
    {
        std::vector<std::vector<int>> localValid;
#pragma omp for schedule(dynamic)
        for (const auto& comb : sumCombinations[key]) {
            if (isValidInRun(run, comb)) {
                localValid.push_back(comb);
            }
        }

#pragma omp critical(combine_valid)
        {
            validCombinations.insert(validCombinations.end(),
                                     localValid.begin(), localValid.end());
        }
    }
    return validCombinations;
}

void KakuroSolver::debugPrintRun(const Run &run) const {
    std::cout << "Run: direction: " << (run.isHorizontal ? "horizontal " : "vertical ") << "sum=" << run.sum
              << ", length=" << run.length << ", cells: ";
    for (const auto &[x, y]: run.cells) {
        std::cout << "(" << x << "," << y << ") ";
    }
    std::cout << "\n";
}

void KakuroSolver::debugPrintAllRuns() const {
    std::cout << "Total runs: " << runs.size() << "\n";
    for (const auto &run: runs) {
        debugPrintRun(run);
    }
}


void KakuroSolver::identifyRuns() {
    runs.clear();
    //std::cout << "\nIdentifying runs...\n";

    // Identify horizontal runs
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (j < size && board[i][j].rightClue > 0) {
                Run run;
                run.sum = board[i][j].rightClue;
                int k = j + 1;
                while (k < size && !board[i][k].isBlack) {
                    run.cells.push_back({i, k});
                    k++;
                }
                run.length = run.cells.size();
                run.isHorizontal = true;
                if (run.length >= 2) {
                    runs.push_back(run);
                }
            }
        }
    }

    // Identify vertical runs
    for (int j = 0; j < size; ++j) {
        for (int i = 0; i < size; ++i) {
            if (i < size && board[i][j].downClue > 0) {
                Run run;
                run.sum = board[i][j].downClue;
                int k = i + 1;
                while (k < size && !board[k][j].isBlack) {
                    run.cells.push_back({k, j});
                    k++;
                }
                run.length = run.cells.size();
                run.isHorizontal = false;
                if (run.length >= 2) {
                    runs.push_back(run);
                }
            }
        }
    }
}

bool KakuroSolver::isValidBoard() const {
    // Track which cells belong to runs
    std::vector<std::vector<std::pair<bool, bool>>> cellInRun(size,
                                                              std::vector<std::pair<bool, bool>>(size, {false,
                                                                                                        false}));  // {horizontal, vertical}
    // check if the board is full of black fields
    bool notBlack = false;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (!board[i][j].isBlack) {
                notBlack = true;
                break;
            }
        }
    }
    if (!notBlack) return false;

    // Check each white cell has both horizontal and vertical clues
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (!board[i][j].isBlack) {
                bool hasHorizontalClue = false;
                bool hasVerticalClue = false;

                // Check horizontal clue
                for (int k = j - 1; k >= 0; k--) {
                    if (board[i][k].isBlack) {
                        hasHorizontalClue = (board[i][k].rightClue > 0);
                        break;
                    }
                }

                // Check vertical clue
                for (int k = i - 1; k >= 0; k--) {
                    if (board[k][j].isBlack) {
                        hasVerticalClue = (board[k][j].downClue > 0);
                        break;
                    }
                }

                if (!hasHorizontalClue || !hasVerticalClue) {
                    return false;
                }
            }
        }
    }
    // Track runs and verify cell participation
    for (const auto &run: runs) {
        for (const auto &[x, y]: run.cells) {
            // Determine if this is a horizontal or vertical run
            bool isHorizontal = true;
            if (run.cells.size() > 1) {
                isHorizontal = (run.cells[0].first == run.cells[1].first);
            }

            if (isHorizontal) {
                cellInRun[x][y].first = true;
            } else {
                cellInRun[x][y].second = true;
            }
        }
    }

    // Verify all white cells participate in both directions
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (!board[i][j].isBlack) {
                if (!cellInRun[i][j].first || !cellInRun[i][j].second) {
                    return false;
                }
            }
        }
    }

    return true;
}

bool KakuroSolver::isSolutionComplete() const {
    bool isComplete = true;
#pragma omp parallel for collapse(2) reduction(&:isComplete)
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (!board[i][j].isBlack && board[i][j].value == 0) {
                isComplete = false;
            }
        }
    }
    return isComplete;
}

bool KakuroSolver::verifyRuns() const {
    bool isValid = true;
#pragma omp parallel for reduction(&:isValid)
    for (const auto& run : runs) {
        int sum = 0;
        std::set<int> used;
        for (const auto& [x, y] : run.cells) {
            if (!board[x][y].isBlack) {
                int value = board[x][y].value;
#pragma omp critical(verify)
                {
                    if (value < 1 || value > 9 || used.count(value)) {
                        isValid = false;
                    }
                    used.insert(value);
                }
                sum += value;
            }
        }
        isValid &= (sum == run.sum);
    }
    return isValid;
}


bool KakuroSolver::isValuePossibleAtPosition(int x, int y, int value) {
    // Check row constraints
    for (int j = 0; j < size; ++j) {
        if (j != y && !board[x][j].isBlack && board[x][j].value == value ||
        j != x && !board[j][y].isBlack && board[j][y].value == value) {
            return false;
        }
    }
    return true;
}

SolveResult KakuroSolver::solve(int runIndex, std::vector<std::vector<Cell>>& solution) {
    if (termination_flag && termination_flag->load()) return SolveResult::NO_SOLUTION;
    if (!isValidRunLengths()) return SolveResult::NO_SOLUTION;

    if (runIndex == runs.size()) {
        if (!isSolutionComplete() || !verifyRuns()) return SolveResult::NO_SOLUTION;

#pragma omp atomic
        solutionCount++;

        if (solutionCount == 1) {
#pragma omp critical
            {
                solution = board;
            }
            return SolveResult::UNIQUE_SOLUTION;
        }
        return SolveResult::MULTIPLE_SOLUTIONS;
    }

    std::vector<std::vector<int>> possibleValues;
#pragma omp critical
    {
        const Run& run = runs[runIndex];
        possibleValues = getPossibleValues(run);
    }

    SolveResult result = SolveResult::NO_SOLUTION;

    if (possibleValues.size() > 8 && runIndex < runs.size() / 2) {
        auto localBoard = board;

#pragma omp parallel for schedule(dynamic) shared(result)
        for (size_t i = 0; i < possibleValues.size(); i++) {
            if (solutionCount > 1) continue;

            std::vector<int> values;
#pragma omp critical
            {
                if (i < possibleValues.size()) {
                    values = possibleValues[i];
                }
            }

            if (values.empty()) continue;

            auto threadBoard = localBoard;
            bool isValid = true;
            std::vector<std::pair<int, int>> modifications;

#pragma omp critical
            {
                const Run& run = runs[runIndex];
                size_t valueIndex = 0;
                for (size_t j = 0; j < run.cells.size() && valueIndex < values.size(); j++) {
                    auto [x, y] = run.cells[j];
                    if (!threadBoard[x][y].isBlack) {
                        int value = values[valueIndex++];
                        if (threadBoard[x][y].value != 0 && threadBoard[x][y].value != value) {
                            isValid = false;
                            break;
                        }
                        if (threadBoard[x][y].value == 0) {
                            threadBoard[x][y].value = value;
                            modifications.push_back({x, y});
                        }
                    }
                }
            }

            if (isValid) {
#pragma omp critical
                {
                    board = threadBoard;
                }

                auto subResult = solve(runIndex + 1, solution);

#pragma omp critical
                {
                    if (subResult == SolveResult::MULTIPLE_SOLUTIONS) {
                        result = SolveResult::MULTIPLE_SOLUTIONS;
                    } else if (subResult == SolveResult::UNIQUE_SOLUTION && result != SolveResult::MULTIPLE_SOLUTIONS) {
                        result = SolveResult::UNIQUE_SOLUTION;
                    }

                    for (auto [x, y] : modifications) {
                        board[x][y].value = 0;
                    }
                }
            }
        }
    } else {
        const Run &run = runs[runIndex];
        for (const auto& values : possibleValues) {

            if (solutionCount > 1) break;
            bool isValid = true;
            std::vector<std::pair<int, int>> modifications;

            size_t valueIndex = 0;
            for (size_t i = 0; i < run.cells.size() && valueIndex < values.size(); i++) {
                auto [x, y] = run.cells[i];
                if (!board[x][y].isBlack) {
                    int value = values[valueIndex++];
                    if (board[x][y].value != 0 && board[x][y].value != value) {
                        isValid = false;
                        break;
                    }
                    if (board[x][y].value == 0) {
                        board[x][y].value = value;
                        modifications.push_back({x, y});
                    }
                }
            }

            if (isValid) {
                auto subResult = solve(runIndex + 1, solution);
                if (subResult == SolveResult::MULTIPLE_SOLUTIONS) {
                    result = SolveResult::MULTIPLE_SOLUTIONS;
                    break;
                }
                if (subResult == SolveResult::UNIQUE_SOLUTION) {
                    result = SolveResult::UNIQUE_SOLUTION;
                }
            }

            for (auto [x, y] : modifications) {
                board[x][y].value = 0;
            }
        }
    }

    return result;
}

bool KakuroSolver::isValidRunLengths() const {
    bool isValid = true;
#pragma omp parallel for reduction(&:isValid)
    for (const Run& run : runs) {
        int whiteCells = 0;
        for (const auto& [x, y] : run.cells) {
            if (!board[x][y].isBlack) whiteCells++;
        }
        isValid &= (whiteCells >= 2);
    }
    return isValid;
}


KakuroSolver::KakuroSolver(int boardSize) : size(boardSize) {
    // Initialize all cells as black by default
    board = std::vector<std::vector<Cell>>(size, std::vector<Cell>(size, Cell(true)));
    precomputeSumCombinations();
}

void KakuroSolver::writeToFile(const std::string &filename, std::vector<std::vector<Cell>> solution) const {
    std::ofstream outFile(filename);
    if (!outFile) {
        throw std::runtime_error("Could not open file for writing: " + filename);
    }

    // Write dimensions
    outFile << size << " " << size << std::endl;

    // Write board
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            const Cell &cell = solution[i][j];
            outFile << std::setw(8) << " ";  // Initial padding
            if (cell.isBlack) {
                if (cell.downClue > 0 && cell.rightClue > 0) {
                    outFile << std::right << cell.downClue << "\\"
                            << std::left << std::setw(2) << cell.rightClue << "     ";
                } else if (cell.downClue > 0) {
                    outFile << std::right << cell.downClue << "\\" << std::left
                            << std::setw(2) << "     ";
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

void KakuroSolver::setCell(int x, int y, bool isBlack, int downClue, int rightClue) {
    if (x >= 0 && x < size && y >= 0 && y < size) {  // Bounds check
        board[x][y] = Cell(isBlack);
        board[x][y].downClue = downClue;
        board[x][y].rightClue = rightClue;
        board[x][y].value = 0;
    }
}

Cell KakuroSolver::getCell(int x, int y) {
    return board[x][y];
}

SolveResult KakuroSolver::solveBoard(std::vector<std::vector<Cell>> &solution) {
    identifyRuns();

    // First validate the board structure before attempting to solve
    if (!isValidBoard()) {
        return SolveResult::INVALID_BOARD;
    }

    // Sort runs by constraint level (fewer possibilities first)
    std::sort(runs.begin(), runs.end(), [this](const Run &a, const Run &b) {
        return getPossibleValues(a).size() < getPossibleValues(b).size();
    });

    return solve(0, solution);
}

void KakuroSolver::printInitialBoard() const {
    std::cout << "\nInitial board configuration:\n";
    printBoard();
}

void KakuroSolver::initializeBoard(const std::vector<std::vector<Cell>> &newBoard) {
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            this->setCell(i, j,
                           newBoard[i][j].isBlack,
                           newBoard[i][j].downClue,
                           newBoard[i][j].rightClue);
        }
    }
}

