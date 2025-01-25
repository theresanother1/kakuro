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
#include <atomic>
#include <mutex>


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
    bool isHorizontal;
};


enum class SolveResult {
    INVALID_BOARD,
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
    std::atomic<bool>* termination_flag = nullptr;
    std::unordered_map<int, std::vector<std::vector < int>>> sumCombinations;


    void precomputeSumCombinations();
    void generateOptimizedCombinations(int targetSum, int length, std::vector <std::vector<int>> &result);

    void printBoard() const;
    bool isValidRunLengths() const;
    bool isValuePossibleAtPosition(int x, int y, int value);
    bool isValidInRun(const Run &run, const std::vector<int> &values);
    std::vector <std::vector<int>> getPossibleValues(const Run &run);
    bool verifyRuns() const;
    void identifyRuns();
    bool isSolutionComplete() const;
    bool isValidBoard() const;
    SolveResult solve(int runIndex, std::vector<std::vector<Cell>> &solution);

    void debugPrintRun(const Run& run) const ;
    void debugPrintAllRuns() const ;

public:
    KakuroSolver(int boardSize);

    void writeToFile(const std::string &filename, std::vector<std::vector<Cell>> solution) const;

    static KakuroSolver readFromFile(const std::string &filename) {
        std::ifstream inFile(filename);
        if (!inFile) {
            throw std::runtime_error("Could not open file for reading: " + filename);
        }
        std::cout << "REad file" << std::endl;

        // Read dimensions
        int rows, cols;
        inFile >> rows;

        std::cout << "Read dimensions "<< rows << " x " << rows << std::endl;
        // Skip the rest of the first line
        std::string dummy;
        std::getline(inFile, dummy);

        KakuroSolver solver(rows);

        std::cout << "REad file, initialized solver " << std::endl;
        // Read board line by line
        for (int i = 0; i < rows; ++i) {
            std::string line;
            std::getline(inFile, line);

            // Skip initial spaces
            size_t pos = 0;
            while (pos < line.length() && line[pos] == ' ') pos++;

            // Process each cell in the line
            for (int j = 0; j < rows; j++) {
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
    Cell getCell(int x, int y);

    SolveResult solveBoard(std::vector <std::vector<Cell>> &solution);

    void printInitialBoard() const;
};


#endif //PART2_KAKUROSOLVER_H
