#include "KakuroSolver.h"

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

std::vector<std::vector<int>> KakuroSolver::getPossibleValues(const Run &run) {
    int key = run.sum * 100 + run.length;
    auto &combinations = sumCombinations[key];
    std::vector<std::vector<int>> validCombinations;

    for (const auto &comb: combinations) {
        if (isValidInRun(run, comb)) {
            validCombinations.push_back(comb);
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
                    /*std::cout << "Found horizontal run at (" << i << "," << j
                              << ") sum=" << run.sum << " length=" << run.length
                              << " cells: ";
                    for (const auto &cell: run.cells) {
                        std::cout << "(" << cell.first << "," << cell.second << ") ";
                    }
                    std::cout << "\n";*/
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
                    /*std::cout << "Found vertical run at (" << i << "," << j
                              << ") sum=" << run.sum << " length=" << run.length
                              << " cells: ";
                    for (const auto &cell: run.cells) {
                        std::cout << "(" << cell.first << "," << cell.second << ") ";
                    }
                    std::cout << "\n";*/
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
                    /*std::cout << "Cell at (" << i << "," << j << " value: " << board[i][j].value
                              << ") missing clue(s). Horizontal: " << hasHorizontalClue
                              << ", Vertical: " << hasVerticalClue << std::endl;*/
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
                    /*std::cout << "Cell at (" << i << "," << j << " value: " << board[i][j].value
                                << ") not part of both horizontal and vertical runs" << std::endl;*/
                    return false;
                }
            }
        }
    }

    return true;
}

bool KakuroSolver::isSolutionComplete() const {
    // Check if all white cells have values
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (!board[i][j].isBlack && board[i][j].value == 0) {
                /*std::cout << "Unfilled cell at (" << i << "," << j << ")" << std::endl;*/
                return false;
            }
        }
    }

    return true;
}

bool KakuroSolver::verifyRuns() const {
    for (const auto &run: runs) {
        int sum = 0;
        std::set<int> used;
        for (const auto &[x, y]: run.cells) {
            if (!board[x][y].isBlack) {
                int value = board[x][y].value;
                if (value < 1 || value > 9 || used.count(value)) {
                    return false;
                }
                used.insert(value);
                sum += value;
            }
        }
        if (sum != run.sum) {
            return false;
        }
    }
    return true;
}

bool KakuroSolver::isValuePossibleAtPosition(int x, int y, int value) {
    // Check row constraints
    for (int j = 0; j < size; j++) {
        if (j != y && !board[x][j].isBlack && board[x][j].value == value) {
            return false;
        }
    }

    // Check column constraints
    for (int i = 0; i < size; i++) {
        if (i != x && !board[i][y].isBlack && board[i][y].value == value) {
            return false;
        }
    }

    return true;
}

SolveResult KakuroSolver::solve(int runIndex, std::vector<std::vector<Cell>> &solution) {
    if (termination_flag && termination_flag->load()) return SolveResult::NO_SOLUTION;

    if (!isValidRunLengths()) {
        //std::cout << " runs not valid " << std::endl;
        return SolveResult::NO_SOLUTION;
    }

    if (runIndex == runs.size()) {
        // Check if all cells are filled and the solution is valid
        if (!isSolutionComplete() || !verifyRuns()) {
            //std::cout << "Solution incomplete/unverified - not all cells filled\n";
            return SolveResult::NO_SOLUTION;
        }

        solutionCount++;
        //std::cout << "\nFound solution #" << solutionCount << ":\n";
        //printBoard();

        if (solutionCount == 1) {
            return SolveResult::UNIQUE_SOLUTION;  // Keep searching to verify uniqueness
        }
        return SolveResult::MULTIPLE_SOLUTIONS;  // Found more than one solution
    }

    const Run &run = runs[runIndex];
    auto possibleValues = getPossibleValues(run);

    // Sort combinations by how many cells they can fill immediately
    std::sort(possibleValues.begin(), possibleValues.end(),
              [this, &run](const std::vector<int>& a, const std::vector<int>& b) {
                  int scoreA = 0, scoreB = 0;
                  for (size_t i = 0; i < run.cells.size(); i++) {
                      const auto& [x, y] = run.cells[i];
                      if (board[x][y].value == 0) {
                          if (isValuePossibleAtPosition(x, y, a[i])) scoreA++;
                          if (isValuePossibleAtPosition(x, y, b[i])) scoreB++;
                      }
                  }
                  return scoreA > scoreB;
              });


    for (const auto &values: possibleValues) {
        bool isValid = true;
        std::vector<std::pair<int, int>> modifications;

        // Try placing values
        size_t valueIndex = 0;
        for (size_t i = 0; i < run.cells.size() && valueIndex < values.size(); i++) {
            auto [x, y] = run.cells[i];
            if (!board[x][y].isBlack) {
                int value = values[valueIndex++];

                // Check if cell already has a different value
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
            auto result = solve(runIndex + 1, solution);
            if (result == SolveResult::MULTIPLE_SOLUTIONS) {
                return SolveResult::MULTIPLE_SOLUTIONS;  // Propagate multiple solutions immediately
            }
            if (result == SolveResult::UNIQUE_SOLUTION && solutionCount > 1) {
                return SolveResult::MULTIPLE_SOLUTIONS;  // Found another solution
            }
        }

        // Backtrack

        for (auto [x, y]: modifications) {
            board[x][y].value = 0;
        }
    }

    if (solutionCount == 1) {
        return SolveResult::UNIQUE_SOLUTION;
    } else if (solutionCount > 1) {
        return SolveResult::MULTIPLE_SOLUTIONS;
    } else {
        //std::cout << "no result at end " << solutionCount << std::endl;
        //printBoard();
        return SolveResult::NO_SOLUTION;
    }
}

bool KakuroSolver::isValidRunLengths() const {
    for (const Run &run: runs) {
        // Count non-black cells in the run
        int whiteCells = 0;
        for (const auto &[x, y]: run.cells) {
            if (!board[x][y].isBlack) {
                whiteCells++;
            }
        }

        // Each run must have at least 2 white cells
        if (whiteCells < 2) {
            return false;
        }
    }
    return true;
}

KakuroSolver::KakuroSolver(int boardSize) : size(boardSize) {
    // Initialize all cells as black by default
    board = std::vector<std::vector<Cell>>(size, std::vector<Cell>(size, Cell(true)));
    precomputeSumCombinations();
    //std::cout << "Board initialized\n";
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

