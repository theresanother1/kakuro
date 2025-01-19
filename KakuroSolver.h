//
// Created by Theresa on 19.01.2025.
//

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


#ifndef PART2_KAKUROSOLVER_H
#define PART2_KAKUROSOLVER_H


struct Cell {
    bool isBlack;
    int value;
    int downClue;
    int rightClue;

    Cell() : isBlack(true), value(0), downClue(0), rightClue(0) {}

    explicit Cell(bool black) : isBlack(black), value(0), downClue(0), rightClue(0) {}
};

struct Run {
    std::vector<std::pair<int, int>> cells;  // Coordinates of cells in this run
    int sum{};
    int length{};
};


enum class SolveResult {
    NO_SOLUTION,
    UNIQUE_SOLUTION,
    MULTIPLE_SOLUTIONS
};


class KakuroSolver {
private:
    std::vector <std::vector<Cell>> board;
    std::vector <Run> runs;
    int size;
    int solutionCount = 0;
    long long backtrackCount = 0;
    long long combinationsTried = 0;
    std::unordered_map<int, std::vector<std::vector < int>>> sumCombinations;


    void precomputeSumCombinations();

    void generateOptimizedCombinations(int targetSum, int length, std::vector <std::vector<int>> &result);

    void printBoard() const;


    bool isValidInRun(const Run &run, const std::vector<int> &values);

    std::vector <std::vector<int>> getPossibleValues(const Run &run);

    /*
    void identifyRuns() {
        runs.clear();
        //std::cout << "\nIdentifying runs...\n";

        // Identify horizontal runs
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                // The bug was here - we were checking rightClue > 0 on the wrong cell!
                // Changed from: if (board[i][j].rightClue > 0) {
                if (j < size && board[i][j].rightClue > 0) {  // Make sure we're in bounds
                    Run run;
                    run.sum = board[i][j].rightClue;
                    int k = j + 1;
                    // Continue collecting cells until we hit a black cell or board edge
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
                        for (const auto& cell : run.cells) {
                            std::cout << "(" << cell.first << "," << cell.second << ") ";
                        }
                        std::cout << "\n";
                    }
                }
            }
        }

        // Identify vertical runs - similar fix needed here
        for (int j = 0; j < size; ++j) {
            for (int i = 0; i < size; ++i) {
                // Changed from: if (board[i][j].downClue > 0) {
                if (i < size && board[i][j].downClue > 0) {  // Make sure we're in bounds
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
                        for (const auto& cell : run.cells) {
                            std::cout << "(" << cell.first << "," << cell.second << ") ";
                        }
                        std::cout << "\n";
                    }
                }
            }
        }

        std::cout << "Total runs found: " << runs.size() << "\n";
    }
*/

    void identifyRuns();

    bool isValidBoard() const;

    /*
    bool isValidBoard() const {
        // Track which cells belong to runs
        std::vector<std::vector<std::pair<bool, bool>>> cellInRun(size,
                                                                  std::vector<std::pair<bool, bool>>(size, {false, false}));  // {horizontal, vertical}


        // check if the board is full of black fields
        bool notBlack = false;
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j< size; ++j) {
                if (!board[i][j].isBlack) {
                    notBlack = true;
                    break;
                }
            }
        }
        if (!notBlack) return false;


        // First identify all runs in the board
        // For each black cell with a clue, mark all white cells it governs
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (board[i][j].isBlack) {
                    // Check right clue
                    if (board[i][j].rightClue > 0) {
                        int k = j + 1;
                        while (k < size && !board[i][k].isBlack) {
                            cellInRun[i][k].first = true;  // Mark horizontal run
                            k++;
                        }
                    }

                    // Check down clue
                    if (board[i][j].downClue > 0) {
                        int k = i + 1;
                        while (k < size && !board[k][j].isBlack) {
                            cellInRun[k][j].second = true;  // Mark vertical run
                            k++;
                        }
                    }
                }
            }
        }

        // Verify all white cells belong to both a horizontal and vertical run
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (!board[i][j].isBlack) {
                    if (!cellInRun[i][j].first || !cellInRun[i][j].second) {
                        std::cout << "Cell at (" << i << "," << j
                                  << ") not part of both a horizontal and vertical run" << std::endl;
                        return false;
                    }
                }
            }
        }

        // Verify all white cells participate in both directions
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (!board[i][j].isBlack) {
                    if (!cellInRun[i][j].first || !cellInRun[i][j].second) {
                        std::cout << "Cell at (" << i << "," << j
                                  << ") not part of both horizontal and vertical runs" << std::endl;
                        return false;
                    }
                }
            }
        }

        return true;
    }
*/
    /*
    bool isValidBoard() const {
        // Track which cells belong to runs
        std::vector<std::vector<std::pair<bool, bool>>> cellInRun(size,
                                                                  std::vector<std::pair<bool, bool>>(size, {false, false}));  // {horizontal, vertical}


        // check if the board is full of black fields
        bool notBlack = false;
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j< size; ++j) {
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
                       std::cout << "Cell at (" << i << "," << j
                                  << ") missing clue(s). Horizontal: " << hasHorizontalClue
                                  << ", Vertical: " << hasVerticalClue << std::endl;
                        return false;
                    }
                }
            }
        }

        // Track runs and verify cell participation
        for (const auto& run : runs) {
            for (const auto& [x, y] : run.cells) {
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
                        std::cout << "Cell at (" << i << "," << j
                                  << ") not part of both horizontal and vertical runs" << std::endl;
                        return false;
                    }
                }
            }
        }

        return true;
    }
*/
    bool isSolutionComplete() const;


    SolveResult solve(int runIndex, std::vector <std::vector<Cell>> &solution);


public:
    KakuroSolver(int boardSize);

    void writeToFile(const std::string &filename, std::vector<std::vector<Cell>> solution) const;
    static KakuroSolver readFromFile(const std::string &filename) {
        std::ifstream inFile(filename);
        if (!inFile) {
            throw std::runtime_error("Could not open file for reading: " + filename);
        }

        // Read dimensions
        int rows, cols;
        inFile >> rows >> cols;
        if (rows != cols) {
            //throw std::runtime_error("Non-square board dimensions not supported");
        }

        // Skip the rest of the first line
        std::string dummy;
        std::getline(inFile, dummy);

        KakuroSolver solver(rows);

        // Read board line by line
        for (int i = 0; i < rows; ++i) {
            std::string line;
            std::getline(inFile, line);

            // Skip initial spaces
            size_t pos = 0;
            while (pos < line.length() && line[pos] == ' ') pos++;

            // Process each cell in the line
            for (int j = 0; j < cols; j++) {
                // Skip spaces between cells
                while (pos < line.length() && line[pos] == ' ') pos++;

                if (pos >= line.length()) break;

                // Parse cell content
                if (line[pos] == '#') {
                    solver.setCell(i, j, true);
                    pos += 8;  // Skip the rest of the cell width
                } else if (line[pos] == '\\') {
                    // Right clue only
                    pos++;  // Skip backslash
                    int rightClue = 0;
                    while (pos < line.length() && std::isdigit(line[pos])) {
                        rightClue = rightClue * 10 + (line[pos] - '0');
                        pos++;
                    }
                    solver.setCell(i, j, true, 0, rightClue);
                } else if (std::isdigit(line[pos])) {
                    int firstNum = 0;
                    while (pos < line.length() && std::isdigit(line[pos])) {
                        firstNum = firstNum * 10 + (line[pos] - '0');
                        pos++;
                    }

                    if (pos < line.length() && line[pos] == '\\') {
                        // Both down and right clues
                        pos++;  // Skip backslash
                        int rightClue = 0;
                        while (pos < line.length() && std::isdigit(line[pos])) {
                            rightClue = rightClue * 10 + (line[pos] - '0');
                            pos++;
                        }
                        solver.setCell(i, j, true, firstNum, rightClue);
                    } else {
                        // Regular number
                        solver.setCell(i, j, false);
                        solver.board[i][j].value = firstNum;
                    }
                } else if (line[pos] == '_') {
                    solver.setCell(i, j, false);
                    pos++;
                }

                // Skip to next cell
                while (pos < line.length() && line[pos] == ' ') pos++;
            }
        }

        inFile.close();
        std::cout << "Board successfully read from " << filename << std::endl;
        std::cout << "Verification - read board state:" << std::endl;
        solver.printBoard();

        return solver;
    }

    void setCell(int x, int y, bool isBlack, int downClue = 0, int rightClue = 0);

    bool solveBoard(std::vector <std::vector<Cell>> &solution);

    void printInitialBoard() const;
};


#endif //PART2_KAKUROSOLVER_H
