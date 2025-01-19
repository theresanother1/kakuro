//
// Created by Theresa on 19.01.2025.
//

#include "KakuroSolver.h"

void KakuroSolver::precomputeSumCombinations() {
    //std::cout << "\nPrecomputing sum combinations...\n";

    // Quick lookup arrays for min/max possible sums
    std::vector<int> minSum(10, 0);  // minSum[length] = minimum possible sum for that length
    std::vector<int> maxSum(10, 0);  // maxSum[length] = maximum possible sum for that length

    for (int len = 1; len <= 9; len++) {
        // Minimum sum: 1+2+3+...+len
        minSum[len] = (len * (len + 1)) / 2;
        // Maximum sum: 9+8+7+...+(9-len+1)
        maxSum[len] = (len * (19 - len)) / 2;
    }

    for (int length = 1; length <= 9; ++length) {
        for (int sum = minSum[length]; sum <= maxSum[length]; ++sum) {
            generateOptimizedCombinations(sum, length, sumCombinations[sum * 100 + length]);
        }
    }

    // Print statistics
    /*for (const auto& [key, combinations] : sumCombinations) {
        int sum = key / 100;
        int length = key % 100;
        std::cout << "Sum " << sum << " Length " << length << ": "
                  << combinations.size() << " combinations\n";
    }*/
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

void KakuroSolver::identifyRuns() {
    runs.clear();
    //std::cout << "\nIdentifying runs...\n";

    // First verify basic board structure
    /* for (int i = 0; i < size; i++) {
         if (!board[i][0].isBlack || !board[0][i].isBlack) {
             throw std::runtime_error("First row and column must be black cells");
         }
     }*/

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
                if (run.length > 0) {
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
                if (run.length > 0) {
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

    //std::cout << "Total runs found: " << runs.size() << "\n";
    std::cout << "Current size: " << size << "\n";

}

bool KakuroSolver::isValidBoard() const {
    // Track which cells belong to runs
    std::vector<std::vector<bool>> cellInRun(size, std::vector<bool>(size, false));

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


    // Identify all runs in the board
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (board[i][j].isBlack) {

                // Check right clue
                if (board[i][j].rightClue > 0) {

                    //Check cell to the right if this is a clue
                    if (j + 1 >= size || board[i][j + 1].isBlack) {
                        /*std::cout << "Invalid right clue at (" << i << "," << j
                                  << "): black cell, board edge or clue immediately after" << std::endl;*/
                        //printBoard();
                        return false;
                    }

                    int k = j + 1;
                    while (k < size && !board[i][k].isBlack) {
                        cellInRun[i][k] = true;
                        k++;
                    }
                }

                // Check down clue
                if (board[i][j].downClue > 0) {

                    //Check cell underneath if this is a clue
                    if (i + 1 >= size || board[i + 1][j].isBlack) {
                        /*std::cout << "Invalid down clue at (" << i << "," << j
                                  << "): black cell, board edge or clue immediately below" << std::endl;*/
                        //printBoard();
                        return false;
                    }

                    int k = i + 1;
                    while (k < size && !board[k][j].isBlack) {
                        cellInRun[k][j] = true;
                        k++;
                    }
                }
            }
        }
    }

    // Verify all white cells belong to at least one run
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (!board[i][j].isBlack && !cellInRun[i][j]) {
                std::cout << "Cell at (" << i << "," << j
                          << ") is not part of any run" << std::endl;
                return false;
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
                std::cout << "Unfilled cell at (" << i << "," << j << ")" << std::endl;
                return false;
            }
        }
    }
    return true;
}

SolveResult KakuroSolver::solve(int runIndex, std::vector<std::vector<Cell>> &solution) {
    combinationsTried++;
    if (combinationsTried % 200000 == 0) {
        std::cout << "Progress: tried " << combinationsTried << " combinations\n";
        //printBoard();
    }

    // First validate the board structure before attempting to solve
    if (runIndex == 0 && !isValidBoard()) {
        //std::cout << "Invalid board structure detected\n";
        return SolveResult::NO_SOLUTION;
    }

    if (runIndex == runs.size()) {
        // Check if all cells are filled and the solution is valid
        if (!isSolutionComplete()) {
            std::cout << "Solution incomplete - not all cells filled\n";
            return SolveResult::NO_SOLUTION;
        }

        // Verify all runs
        for (const Run &run: runs) {
            int sum = 0;
            for (const auto &[x, y]: run.cells) {
                if (!board[x][y].isBlack) {
                    sum += board[x][y].value;
                }
            }
            if (sum != run.sum) return SolveResult::NO_SOLUTION;
        }

        solutionCount++;
        std::cout << "\nFound solution #" << solutionCount << ":\n";
        printBoard();

        if (solutionCount == 1) {
            solution = board;
            return SolveResult::UNIQUE_SOLUTION;  // Keep searching to verify uniqueness
        }
        return SolveResult::MULTIPLE_SOLUTIONS;  // Found more than one solution
    }

    const Run &run = runs[runIndex];
    auto possibleValues = getPossibleValues(run);

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

                // Check row and column conflicts
                for (int k = 0; k < size; k++) {
                    if (k != y && !board[x][k].isBlack && board[x][k].value == value) {
                        isValid = false;
                        break;
                    }
                    if (k != x && !board[k][y].isBlack && board[k][y].value == value) {
                        isValid = false;
                        break;
                    }
                }

                if (!isValid) break;

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
        backtrackCount++;
        for (auto [x, y]: modifications) {
            board[x][y].value = 0;
        }
    }

    return solutionCount == 1 ? SolveResult::UNIQUE_SOLUTION : SolveResult::NO_SOLUTION;
}

KakuroSolver::KakuroSolver(int boardSize) : size(boardSize) {
    //std::cout << "Initializing " << boardSize << "x" << boardSize << " board\n";
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
        board[x][y].isBlack = isBlack;
        board[x][y].downClue = downClue;
        board[x][y].rightClue = rightClue;
        board[x][y].value = 0;
    }
}

bool KakuroSolver::solveBoard(std::vector<std::vector<Cell>> &solution) {
    identifyRuns();

    // Sort runs by constraint level (fewer possibilities first)
    std::sort(runs.begin(), runs.end(), [this](const Run &a, const Run &b) {
        return getPossibleValues(a).size() < getPossibleValues(b).size();
    });

    /*std::cout << "\nRuns sorted by constraint level:\n";
    for (const auto& run : runs) {
        std::cout << "Run sum=" << run.sum << " length=" << run.length
                  << " possibilities=" << getPossibleValues(run).size() << "\n";
    }*/

    solutionCount = 0;
    backtrackCount = 0;
    combinationsTried = 0;

    auto result = solve(0, solution);

    switch (result) {
        case SolveResult::NO_SOLUTION:
            //std::cout << "No solution exists\n";
            break;
        case SolveResult::UNIQUE_SOLUTION:
            //std::cout << "Found unique solution!\n";
            // solution vector contains the unique solution
            break;
        case SolveResult::MULTIPLE_SOLUTIONS:
            //std::cout << "Puzzle has multiple solutions - not valid!\n";
            break;
    }

    std::cout //<< "\nSolving completed:\n"
            << "Solutions found: " << solutionCount << "\n" << std::endl;
    //<< "Combinations tried: " << combinationsTried << "\n"
    //<< "Backtracks: " << backtrackCount << "\n";

    return result == SolveResult::UNIQUE_SOLUTION;
}

void KakuroSolver::printInitialBoard() const {
    std::cout << "\nInitial board configuration:\n";
    printBoard();
}
