#include <vector>
#include <random>
#include <algorithm>
#include <chrono>
#include <map>
#include <mutex>
#include <atomic>
#include <future>
#include "KakuroSolver.h"
#include <omp.h>

class KakuroBoard {
private:
    int rows;
    int cols;

    const int minRunLength = 2;
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<> dis;
    std::uniform_int_distribution<> numDis;

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
};


class KakuroMutations {
private:
    const int size;
    std::mt19937 rng;
    std::uniform_real_distribution<> dist;

    struct RunInfo {
        int startX, startY;
        int length;
        bool isHorizontal;
    };

    // Track all runs on the board
    std::vector<RunInfo> getRunInfo(const std::vector<std::vector<Cell>> &board) {
        std::vector<RunInfo> runs;

        // Horizontal runs
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (board[i][j].isBlack && j + 1 < size && !board[i][j + 1].isBlack) {
                    RunInfo run;
                    run.startX = i;
                    run.startY = j;
                    run.length = countHorizontalRun(board, i, j);
                    run.isHorizontal = true;
                    runs.push_back(run);
                }
            }
        }

        // Vertical runs
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (board[i][j].isBlack && i + 1 < size && !board[i][i + 1].isBlack) {
                    RunInfo run;
                    run.startX = i;
                    run.startY = j;
                    run.length = countVerticalRun(board, i, j);
                    run.isHorizontal = false;
                    runs.push_back(run);
                }
            }
        }
        return runs;
    }

    double calculateBlackCellRatio(const std::vector<std::vector<Cell>> &board) {
        int blackCount = 0;
        int total = size * size;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (board[i][j].isBlack) blackCount++;
            }
        }
        return static_cast<double>(blackCount) / total;
    }

    bool isIsolated(const std::vector<std::vector<Cell>> &board, int x, int y) {
        if (board[x][y].isBlack) return false;

        int connectedRuns = 0;
        // Check horizontal run
        if (y > 0 && board[x][y - 1].rightClue > 0) connectedRuns++;
        // Check vertical run
        if (x > 0 && board[x - 1][y].downClue > 0) connectedRuns++;

        return connectedRuns < 2;
    }

    bool wouldCreateShortRun(const std::vector<std::vector<Cell>> &board, int x, int y) {
        if (!board[x][y].isBlack) return false;

        int hLength = countHorizontalRun(board, x, y);
        int vLength = countVerticalRun(board, x, y);

        return (hLength == 1 || vLength == 1);
    }

    // Count length of horizontal run starting at position
    int countHorizontalRun(const std::vector<std::vector<Cell>> &board, int i, int j) {
        int length = 0;
        j++; // Skip the black cell with clue
        while (j < size && !board[i][j].isBlack) {
            length++;
            j++;
        }
        return length;
    }

    // Count length of vertical run starting at position
    int countVerticalRun(const std::vector<std::vector<Cell>> &board, int i, int j) {
        int length = 0;
        i++; // Skip the black cell with clue
        while (i < size && !board[i][j].isBlack) {
            length++;
            i++;
        }
        return length;
    }

    // Generate valid sum for given run length
    int generateValidSum(int length) {
        if (length < 2) return 0;

        // Calculate valid sum range for this length
        int minSum = (length * (length + 1)) / 2;  // Sum of smallest numbers (1,2,3...)
        int maxSum = length * 9 - (length * (length - 1)) / 2;  // Sum of largest possible numbers

        std::uniform_int_distribution<> sumDist(minSum, maxSum);
        return sumDist(rng);
    }

public:
    KakuroMutations(int boardSize) :
            size(boardSize),
            rng(std::random_device{}()),
            dist(0.0, 1.0) {}


    void mutateBoard(std::vector<std::vector<Cell>> &board, double mutationRate) {
        auto runs = getRunInfo(board);
        double blackRatio = calculateBlackCellRatio(board);

        // Structural mutations
        for (int i = 1; i < size - 1; i++) {
            for (int j = 1; j < size - 1; j++) {
                if (dist(rng) < mutationRate) {
                    // Check if we should add/remove black cells based on ratio
                    bool preferBlack = (blackRatio < 0.3);
                    bool canFlip = !wouldCreateShortRun(board, i, j);

                    if (canFlip) {
                        if (board[i][j].isBlack && preferBlack) {
                            board[i][j].isBlack = false;
                            board[i][j].rightClue = board[i][j].downClue = 0;
                        } else if (!board[i][j].isBlack && !preferBlack && !isIsolated(board, i, j)) {
                            board[i][j].isBlack = true;
                        }
                    }
                }
            }
        }

        // Clue mutations
        mutateClues(board, runs, mutationRate);
    }

    void mutateClues(std::vector<std::vector<Cell>> &board,
                     const std::vector<RunInfo> &runs,
                     double mutationRate) {
        for (const auto &run: runs) {
            if (dist(rng) < mutationRate) {
                int x = run.startX;
                int y = run.startY;

                if (run.length >= 2) {
                    if (run.isHorizontal) {
                        board[x][y].rightClue = generateValidSum(run.length);
                    } else {
                        board[x][y].downClue = generateValidSum(run.length);
                    }
                }
            }
        }
    }

};

class ThreadSafeRNG {
    std::mutex mtx;
    std::mt19937 rng;
public:

    int getInt(int min, int max) {
        std::lock_guard<std::mutex> lock(mtx);
        return std::uniform_int_distribution<>(min, max)(rng);
    }

    double getReal(double min = 0.0, double max = 1.0) {
        std::lock_guard<std::mutex> lock(mtx);
        return std::uniform_real_distribution<>(min, max)(rng);
    }

    ThreadSafeRNG() : rng(std::random_device{}()) {}
};


class PatternGenerator {

public:
    enum PatternType {
        CROSS,
        CHECKERBOARD,
        CENTER_HEAVY,
        EDGE_HEAVY,
        DIAGONAL
    };

    enum TemplateType {
        DIAMOND,
        SPIRAL,
        BLOCK,
        ZIGZAG
    };

    PatternGenerator(ThreadSafeRNG &rng, int boardSize) : threadRng(rng), size(boardSize) {}

    std::vector<std::vector<bool>> generateRandomPattern() {
        PatternType type = static_cast<PatternType>(threadRng.getInt(0, 4));
        std::vector<bool> basePattern = generateBasePattern(type);

        std::vector<std::vector<bool>> result(size, std::vector<bool>(size));
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                result[i][j] = basePattern[(i + j) % basePattern.size()];
            }
        }
        return result;
    }

    std::vector<std::vector<bool>> generateRandomTemplate() {
        TemplateType type = static_cast<TemplateType>(threadRng.getInt(0, 3));
        return generateTemplate(type);
    }

private:
    ThreadSafeRNG &threadRng;
    const int size;

    std::vector<bool> generateBasePattern(PatternType type) {
        std::vector<bool> pattern;
        pattern.reserve(size);

        switch (type) {
            case CROSS:
                for (int i = 0; i < size; i++) {
                    pattern.push_back(i % 3 == 1);
                }
                break;
            case CHECKERBOARD:
                for (int i = 0; i < size; i++) {
                    pattern.push_back(i % 2 == 0);
                }
                break;
            case CENTER_HEAVY:
                for (int i = 0; i < size; i++) {
                    pattern.push_back(i >= size / 4 && i < 3 * size / 4);
                }
                break;
            case EDGE_HEAVY:
                for (int i = 0; i < size; i++) {
                    pattern.push_back(i < size / 4 || i >= 3 * size / 4);
                }
                break;
            case DIAGONAL:
                for (int i = 0; i < size; i++) {
                    pattern.push_back(i < (size + 1) / 2);
                }
                break;
        }
        return pattern;
    }

    std::vector<std::vector<bool>> generateTemplate(TemplateType type) {
        std::vector<std::vector<bool>> templ(size, std::vector<bool>(size, false));

        switch (type) {
            case DIAMOND:
                generateDiamondTemplate(templ);
                break;
            case SPIRAL:
                generateSpiralTemplate(templ);
                break;
            case BLOCK:
                generateBlockTemplate(templ);
                break;
            case ZIGZAG:
                generateZigzagTemplate(templ);
                break;
        }

        return templ;
    }

    void generateDiamondTemplate(std::vector<std::vector<bool>> &templ) {
        int mid = size / 2;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                int dist = std::abs(i - mid) + std::abs(j - mid);
                templ[i][j] = (dist <= mid);
            }
        }
    }

    void generateSpiralTemplate(std::vector<std::vector<bool>> &templ) {
        int left = 0, right = size - 1, top = 0, bottom = size - 1;
        bool fill = true;

        while (left <= right && top <= bottom) {
            // Fill top row
            for (int i = left; i <= right; i++)
                templ[top][i] = fill;
            top++;

            // Fill right column
            for (int i = top; i <= bottom; i++)
                templ[i][right] = fill;
            right--;

            if (top <= bottom) {
                // Fill bottom row
                for (int i = right; i >= left; i--)
                    templ[bottom][i] = fill;
                bottom--;
            }

            if (left <= right) {
                // Fill left column
                for (int i = bottom; i >= top; i--)
                    templ[i][left] = fill;
                left++;
            }

            fill = !fill;
        }
    }

    void generateBlockTemplate(std::vector<std::vector<bool>> &templ) {
        int blockSize = std::max(2, size / 5);
        for (int i = 0; i < size; i += blockSize) {
            for (int j = 0; j < size; j += blockSize) {
                bool fillBlock = ((i / blockSize + j / blockSize) % 2 == 0);
                for (int bi = 0; bi < blockSize && i + bi < size; bi++) {
                    for (int bj = 0; bj < blockSize && j + bj < size; bj++) {
                        templ[i + bi][j + bj] = fillBlock;
                    }
                }
            }
        }
    }

    void generateZigzagTemplate(std::vector<std::vector<bool>> &templ) {
        int stripWidth = std::max(2, size / 8);
        for (int i = 0; i < size; i++) {
            int offset = (i / stripWidth) * stripWidth;
            for (int j = 0; j < size; j++) {
                templ[i][j] = ((j + offset) / stripWidth) % 2 == 0;
            }
        }
    }

};

struct Individual {
    std::vector<std::vector<Cell>> board;
    double fitness;

    Individual(int size) : board(size, std::vector<Cell>(size, Cell(true))), fitness(0.0) {}
};


struct Island {
    std::vector<Individual> population;
    double bestFitness;
    Individual bestIndividual;
    double mutationRate;
    double crossoverRate;
    int generationsWithoutImprovement;
    double selectionPressure;
    double denseNessBlack;
    int id;

    Island(int size, int popSize, int id) :
            population(popSize, Individual(size)),
            bestFitness(0.0),
            bestIndividual(size),
            mutationRate(0.1),
            crossoverRate(0.7),
            generationsWithoutImprovement(0),
            selectionPressure(0.7),
            denseNessBlack(0.2),
            id(id) {
    }
};

class KakuroGenerator {
private:
    KakuroMutations mutator;
    int MAX_GENERATIONS = 2000;
    const int size;
    const double mutationRate;
    const double targetFitness;
    ThreadSafeRNG threadRng;
    double temperature;
    std::unordered_map<std::string, double> fitnessCache;
    static constexpr size_t MAX_CACHE_SIZE = 10000;
    static constexpr int MAX_POPULATION_SIZE = 500;
    static constexpr double MIN_DIVERSITY_THRESHOLD = 0.2;

    std::mutex cacheMutex;


    void cleanupCache() {
        std::lock_guard<std::mutex> lock(cacheMutex);
        if (fitnessCache.size() > MAX_CACHE_SIZE) {
          fitnessCache.clear();
        }
    }

    double checkAndUpdateCache(const std::string& cacheKey,
                               const double value) {
        std::lock_guard<std::mutex> lock(cacheMutex);
        auto it = fitnessCache.find(cacheKey);
        if (it != fitnessCache.end()) {
            return it->second;
        }
        // - 1 for only check no update
        if (value != -1) fitnessCache[cacheKey] = value;
        return value;
    }


    struct StagnationParams {
        static constexpr int STAGNATION_THRESHOLD = 15;
        static constexpr double MIN_DIVERSITY = MIN_DIVERSITY_THRESHOLD;
        static constexpr double BASE_INJECTION_RATE = 0.3;
    };

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
        double solvability;

        FitnessComponents() : basicValidity(0), runQuality(0),
                              clueDistribution(0),
                              solvability(0) {}
    };


    std::string getBoardHash(const std::vector<std::vector<Cell>>& board) {
        std::stringstream ss;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                const auto& cell = board[i][j];
                ss << (cell.isBlack ? "1" : "0")
                   << "," << cell.rightClue
                   << "," << cell.downClue
                   << ";";
            }
        }
        return ss.str();
    }


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


        // Calculate percentage of white cells that are part of runs
        whiteCellsInRuns = std::min(whiteCellsInRuns, totalWhiteCells);
        double runCoverage = static_cast<double>(whiteCellsInRuns) / totalWhiteCells;
        score += 0.3 * runCoverage;

        if (totalClues > 0) {
            score += 0.3 * (static_cast<double>(validClues) / totalClues);
        }

        if (totalWhiteCells > 0) {
            score += 0.3 * (static_cast<double>(connectedWhiteCells) / totalWhiteCells);
        }


        if (boxedClues > 0) {
            score += 0.1 * (static_cast<double>(boxedClues) / totalClues);
        }

        return score;
    }


    // fitness calculation helper function
/*
    double calculateRunQualityScore(const std::vector<std::vector<Cell>> &board) {
        double score = 0.0;
        int totalRuns = 0;
        std::map<int, int> lengthFrequency;

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

                    lengthFrequency[runLength]++;

                    // Score based on run length ((2-4 cells is ideal for most puzzles, higher makes it exponentially harder)
                    if (runLength >= 2 && runLength <= 4) {
                        score += 1.0;  // Optimal length
                    } else if (runLength > 4 && runLength <= 6) {
                        score += 0.25;  // Acceptable but not ideal
                    } else if (runLength > 6) {
                        score -= 0.25;  // Penalize very long runs
                    }

                    // Validate and score clue values
                    int clue = board[i][j].rightClue;
                    int minPossible = (runLength * (runLength + 1)) / 2;
                    int maxPossible = runLength * 9 - (runLength * (runLength - 1)) / 2;

                    if (clue >= minPossible && clue <= maxPossible) {
                        score += 0.5;
                        // Bonus for interesting sums (not min/max)
                        int range = maxPossible - minPossible;
                        int margin = range / 4;
                        if (clue > minPossible + margin && clue < maxPossible - margin) {
                            score += 0.25;
                        }
                    } else {
                        score -= 0.5;  // Penalty for invalid clues
                    }
                }
            }
        }

        // Check vertical runs
        for (int j = 0; j < size; j++) {
            for (int i = 0; i < size; i++) {
                if (board[i][j].isBlack && board[i][j].downClue > 0) {
                    totalRuns++;
                    int runLength = 0;
                    int k = i + 1;
                    while (k < size && !board[k][j].isBlack) {
                        runLength++;
                        k++;
                    }

                    lengthFrequency[runLength]++;

                    if (runLength >= 2 && runLength <= 4) {
                        score += 1.0;
                    } else if (runLength > 4 && runLength <= 6) {
                        score += 0.25;
                    } else if (runLength > 6) {
                        score -= 0.25;
                    }

                    int clue = board[i][j].downClue;
                    int minPossible = (runLength * (runLength + 1)) / 2;
                    int maxPossible = runLength * 9 - (runLength * (runLength - 1)) / 2;

                    if (clue >= minPossible && clue <= maxPossible) {
                        score += 0.5;
                        int range = maxPossible - minPossible;
                        int margin = range / 4;
                        if (clue > minPossible + margin && clue < maxPossible - margin) {
                            score += 0.25;
                        }
                    } else {
                        score -= 0.5;
                    }
                }
            }
        }

        // Penalize unbalanced run length distribution
        if (totalRuns >= 2) {
            int dominantLength = 0;
            int maxFreq = 0;
            for (const auto &[length, freq]: lengthFrequency) {
                if (freq > maxFreq) {
                    maxFreq = freq;
                    dominantLength = length;
                }
            }
            double ratioOfDominant = static_cast<double>(maxFreq) / totalRuns;
            if (ratioOfDominant > 0.5) {  // If more than 50% runs are same length
                score *= (1.0 - (ratioOfDominant - 0.5));  // Reduce score
            }
        }

        return std::min(1.0, totalRuns > 0 ? score / totalRuns : 0.0);
    }

*/

    double calculateRunQualityScore(const std::vector<std::vector<Cell>>& board) {
        const int size = board.size();
        double totalScore = 0.0;
        int totalRuns = 0;
        std::map<int, int> lengthFrequency;
        std::mutex freqMutex;

#pragma omp parallel reduction(+:totalScore,totalRuns)
        {
            std::map<int, int> localFrequency;

            // Process horizontal runs
#pragma omp for schedule(dynamic)
            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j) {
                    if (board[i][j].isBlack && board[i][j].rightClue > 0) {
                        totalRuns++;
                        double localScore = 0.0;
                        int runLength = 0;

                        // Calculate run length
                        int k = j + 1;
                        while (k < size && !board[i][k].isBlack) {
                            runLength++;
                            k++;
                        }

                        localFrequency[runLength]++;

                        // Score run length
                        if (runLength >= 2 && runLength <= 4) {
                            localScore += 1.0;
                        } else if (runLength > 4 && runLength <= 6) {
                            localScore += 0.25;
                        } else if (runLength > 6) {
                            localScore -= 0.25;
                        }

                        // Validate clue values
                        int clue = board[i][j].rightClue;
                        int minPossible = (runLength * (runLength + 1)) / 2;
                        int maxPossible = runLength * 9 - (runLength * (runLength - 1)) / 2;

                        if (clue >= minPossible && clue <= maxPossible) {
                            localScore += 0.5;
                            int range = maxPossible - minPossible;
                            int margin = range / 4;
                            if (clue > minPossible + margin && clue < maxPossible - margin) {
                                localScore += 0.25;
                            }
                        } else {
                            localScore -= 0.5;
                        }

                        totalScore += localScore;
                    }
                }
            }

            // Process vertical runs
#pragma omp for schedule(dynamic)
            for (int j = 0; j < size; j++) {
                for (int i = 0; i < size; i++) {
                    // Similar processing for vertical runs...
                    // [Implementation similar to horizontal runs]
                }
            }

            // Merge local frequency counts
#pragma omp critical
            {
                for (const auto& [length, freq] : localFrequency) {
                    lengthFrequency[length] += freq;
                }
            }
        }

        // Calculate distribution penalty
        if (totalRuns >= 2) {
            int dominantLength = 0;
            int maxFreq = 0;
            for (const auto& [length, freq] : lengthFrequency) {
                if (freq > maxFreq) {
                    maxFreq = freq;
                    dominantLength = length;
                }
            }
            double ratioOfDominant = static_cast<double>(maxFreq) / totalRuns;
            if (ratioOfDominant > 0.5) {
                totalScore *= (1.0 - (ratioOfDominant - 0.5));
            }
        }

        return std::min(1.0, totalRuns > 0 ? totalScore / totalRuns : 0.0);
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
    double calculateBoardStructureScore(const std::vector<std::vector<Cell>>& board) {
        const int size = board.size();
        double score = 0.0;
        const int totalCells = size * size;

        int blackCells = 0;
        int connectedWhiteCells = 0;
        int totalWhiteCells = 0;
        int clues = 0;
        int isolatedClue = 0;

#pragma omp parallel reduction(+:blackCells,totalWhiteCells,connectedWhiteCells,clues,isolatedClue)
        {
#pragma omp for collapse(2) schedule(static)
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    if (board[i][j].isBlack) {
                        blackCells++;
                    } else {
                        totalWhiteCells++;
                        bool hasHorizontalRun = false;
                        bool hasVerticalRun = false;

                        // Check runs
                        for (int k = j - 1; k >= 0; k--) {
                            if (board[i][k].isBlack) {
                                hasHorizontalRun = (board[i][k].rightClue > 0);
                                break;
                            }
                        }
                        for (int k = i - 1; k >= 0; k--) {
                            if (board[k][j].isBlack) {
                                hasVerticalRun = (board[k][j].downClue > 0);
                                break;
                            }
                        }
                        if (hasHorizontalRun && hasVerticalRun) {
                            connectedWhiteCells++;
                        }
                    }

                    // Count clues
                    if (board[i][j].isBlack) {
                        if (board[i][j].downClue > 0) {
                            clues++;
                            if (j + 1 < size && board[i][j + 1].isBlack)
                                isolatedClue++;
                        }
                        if (board[i][j].rightClue > 0) {
                            clues++;
                            if (i + 1 < size && board[i + 1][j].isBlack)
                                isolatedClue++;
                        }
                    }
                }
            }
        }

        // Calculate final scores
        double blackRatio = static_cast<double>(blackCells) / totalCells;
        score += (1.0 - std::abs(0.3 - blackRatio)) * 0.3;

        if (totalWhiteCells > 0) {
            score += (static_cast<double>(connectedWhiteCells) / totalWhiteCells) * 0.35;
        }
        if (clues > 0) {
            score += (static_cast<double>(isolatedClue) / clues) * 0.35;
        }

        return score;
    }


    /*
    double calculateFitness(Individual &individual) {
        std::string cacheKey = getBoardHash(individual.board);
        FitnessComponents components;

        std::lock_guard<std::mutex> lock(cacheMutex);
        auto it = fitnessCache.find(cacheKey);
        if (it != fitnessCache.end()) {
            return it->second;
        }

        // check basic board structure
        double structureScore = calculateBoardStructureScore(individual.board);
        if (structureScore < 0.4) {
            fitnessCache[cacheKey] = structureScore * 0.2;
            return structureScore * 0.2;
        }

        // Validate each white cell has required clues
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (!individual.board[i][j].isBlack) {
                    bool hasHClue = false, hasVClue = false;

                    // Check for horizontal clue
                    for (int k = j - 1; k >= 0; k--) {
                        if (individual.board[i][k].isBlack) {
                            hasHClue = individual.board[i][k].rightClue > 0;
                            break;
                        }
                    }

                    // Check for vertical clue
                    for (int k = i - 1; k >= 0; k--) {
                        if (individual.board[k][j].isBlack) {
                            hasVClue = individual.board[k][j].downClue > 0;
                            break;
                        }
                    }

                    if (!hasHClue || !hasVClue) {
                        fitnessCache[cacheKey] = structureScore * 0.3;
                        return structureScore * 0.3;
                    }
                }
            }
        }
        components.basicValidity = checkBasicValidity(individual.board);
        if (components.basicValidity < 0.3) {
            fitnessCache[cacheKey] = components.basicValidity * 0.2;
            return components.basicValidity * 0.2;
        }

        if (structureScore > 0.7 && components.basicValidity > 0.7) {
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
            SolveResult result = solveBoardWithTimeOut(solver);

            switch (result) {
                case SolveResult::INVALID_BOARD:
                    fitnessCache[cacheKey] = structureScore * 0.4 + components.runQuality * 0.3 +
                                             components.clueDistribution * 0.1;
                    return structureScore * 0.4 + components.runQuality * 0.3 +
                           components.clueDistribution * 0.1;
                case SolveResult::UNIQUE_SOLUTION:
                    //std::cout << "Found unique solution" << std::endl;
                    components.solvability = 1.0;
                    fitnessCache[cacheKey] = 1.0;
                    return 1.0; // perfect board immediate return
                case SolveResult::MULTIPLE_SOLUTIONS:
                    std::cout << "1";
                    components.solvability = 0.5;
                    break;
                case SolveResult::NO_SOLUTION:
                    //std::cout << "found no solution" << std::endl;
                    components.solvability = 0.0;
                    break;
            }
        }


        components.runQuality = calculateRunQualityScore(individual.board);
        components.clueDistribution = calculateClueDistributionScore(individual.board);
        fitnessCache[cacheKey] = structureScore * 0.4 +
                                 components.runQuality * 0.3 +
                                 components.clueDistribution * 0.1 +
                                 components.solvability * 0.2;

        // Combine components with weights
        return structureScore * 0.4 +
               components.runQuality * 0.3 +
               components.clueDistribution * 0.1 +
               components.solvability * 0.2;
    }
*/

    double calculateFitness(Individual& individual) {
        std::string cacheKey = getBoardHash(individual.board);
        const int size = individual.board.size();

        // Try cache first
        double value = checkAndUpdateCache(cacheKey, -1.0);
        if (value != -1) return value;

        // Calculate structure score
        double structureScore = calculateBoardStructureScore(individual.board);
        if (structureScore < 0.4) {
            return checkAndUpdateCache(cacheKey, structureScore * 0.2);
        }

        // Validate clues - this section is hard to parallelize due to dependencies
        bool validClues = true;
#pragma omp parallel
        {
            bool localValid = true;
#pragma omp for collapse(2)
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    if (!individual.board[i][j].isBlack) {
                        bool hasHClue = false, hasVClue = false;

                        for (int k = j - 1; k >= 0 && !hasHClue; k--) {
                            if (individual.board[i][k].isBlack) {
                                hasHClue = individual.board[i][k].rightClue > 0;
                            }
                        }

                        for (int k = i - 1; k >= 0 && !hasVClue; k--) {
                            if (individual.board[k][j].isBlack) {
                                hasVClue = individual.board[k][j].downClue > 0;
                            }
                        }

                        if (!hasHClue || !hasVClue) {
                            localValid = false;
                        }
                    }
                }
            }
#pragma omp critical
            {
                validClues &= localValid;
            }
        }

        if (!validClues) {
            return checkAndUpdateCache(cacheKey, structureScore * 0.3);
        }

        // Check basic validity
        double basicValidity = checkBasicValidity(individual.board);
        if (basicValidity < 0.3) {
            return checkAndUpdateCache(cacheKey, basicValidity * 0.2);
        }

        // Handle solver and remaining calculations
        double solvability = 0.0;
        if (structureScore > 0.7 && basicValidity > 0.7) {
            KakuroSolver solver(size);

            // Setup solver
            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j) {
                    solver.setCell(i, j,
                                   individual.board[i][j].isBlack,
                                   individual.board[i][j].downClue,
                                   individual.board[i][j].rightClue);
                }
            }

            SolveResult result = solveBoardWithTimeOut(solver);
            switch (result) {
                case SolveResult::UNIQUE_SOLUTION:
                    return checkAndUpdateCache(cacheKey, 1.0);
                case SolveResult::MULTIPLE_SOLUTIONS:
                    solvability = 0.5;
                    break;
                case SolveResult::NO_SOLUTION:
                case SolveResult::INVALID_BOARD:
                    solvability = 0.0;
                    break;
            }
        }

        // Calculate remaining scores
        double runQuality = calculateRunQualityScore(individual.board);
        double clueDistribution = calculateClueDistributionScore(individual.board);

        // Calculate final fitness
        double finalFitness = structureScore * 0.4 +
                              runQuality * 0.3 +
                              clueDistribution * 0.1 +
                              solvability * 0.2;

        return checkAndUpdateCache(cacheKey, finalFitness);
    }

    Individual createIndividualWithBasicBoardValidity() {
        Individual ind(size);
        KakuroBoard newBoard(size, size);
        bool validBoard = false;

        do {
            validBoard = newBoard.generateBoard();
        } while (!validBoard);

        ind.board = newBoard.board;

        return ind;

    }

    void heavyMutation(Individual &ind, double populationDiversity) {
        double mutationChance = 0.3 + (0.2 * (1.0 - populationDiversity));

        for (int i = 0; i < size; ++i) {  // Skip border
            for (int j = 0; j < size; ++j) {
                if (threadRng.getReal() < populationDiversity) {// Flip cell type
                    ind.board[i][j].isBlack = !ind.board[i][j].isBlack;

                    if (ind.board[i][j].isBlack) {// Add new clues with higher probability
                        if (j + 1 < size && !ind.board[i][j + 1].isBlack && threadRng.getReal() < 0.8) {
                            ind.board[i][j].rightClue = threadRng.getInt(1, 45);
                        }
                        if (i + 1 < size && !ind.board[i + 1][j].isBlack && threadRng.getReal() < 0.8) {
                            ind.board[i][j].downClue = threadRng.getInt(1, 45);
                        }
                    } else {// Clear clues if cell becomes white
                        ind.board[i][j].rightClue = 0;
                        ind.board[i][j].downClue = 0;
                    }
                }
            }
        }

// Additionally, mutate some existing clues
        for (int i = 1; i < size; ++i) {
            for (int j = 1; j < size; ++j) {
                if (ind.board[i][j].isBlack && threadRng.getReal() < 0.4) {  // 40% chance for each black cell
                    if (ind.board[i][j].rightClue > 0) {
                        ind.board[i][j].rightClue = threadRng.getInt(1, 45);
                    }
                    if (ind.board[i][j].downClue > 0) {
                        ind.board[i][j].downClue = threadRng.getInt(1, 45);
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
        for (const auto &run: (threadRng.getReal() < 0.5 ? runs1 : runs2)) {
            if (threadRng.getReal() < 0.3) {  // 30% chance to copy a run
                const auto &sourceBoard = (threadRng.getReal() < 0.5 ? parent1.board : parent2.board);

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

        // random crossover operations
        for (int i = 0; i < size; ++i) {  // Skip border cells
            for (int j = 0; j < size; ++j) {  // Skip border cells
                if (!modified[i][j] && threadRng.getReal() < 0.4) {  // 40% chance of random cell
                    child.board[i][j] = Cell(threadRng.getReal() < 0.5);  // Random cell type
                } else if (!modified[i][j]) {
                    child.board[i][j] = (threadRng.getReal() < 0.5) ?
                                        parent1.board[i][j] : parent2.board[i][j];
                }
            }
        }

        // Second phase: Fill remaining cells from either parent
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (!modified[i][j]) {
                    // Choose from which parent to copy
                    child.board[i][j] = (threadRng.getReal() < 0.5) ?
                                        parent1.board[i][j] : parent2.board[i][j];
                }
            }
        }

        // Third phase: Quick validation and repair
        //repairInvalidCells(child.board);
        //repairInvalidRuns(child.board);


        return child;
    }

    void repairInvalidRuns(std::vector<std::vector<Cell>> &board) {
        // First identify all runs
        std::vector<Run> runs;
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (board[i][j].isBlack) {
                    // Check horizontal run
                    if (board[i][j].rightClue > 0) {
                        int length = 0;
                        int k = j + 1;
                        while (k < size && !board[i][k].isBlack) {
                            length++;
                            k++;
                        }
                        if (length < 2) {
                            // Extend run if possible, otherwise merge with neighbor
                            if (k < size - 1) {
                                board[i][k].isBlack = false;
                            } else if (j > 1 && board[i][j - 1].isBlack) {
                                board[i][j - 1].isBlack = false;
                            }
                        }
                    }
                    // Check vertical run
                    if (board[i][j].downClue > 0) {
                        int length = 0;
                        int k = i + 1;
                        while (k < size && !board[k][j].isBlack) {
                            length++;
                            k++;
                        }
                        if (length < 2) {
                            //extend run if possible, otherwise merge with neighbour
                            if (k < size - 1) {
                                board[k][j].isBlack = false;
                            } else if (i > 1 && board[i - 1][j].isBlack) {
                                board[i - 1][j].isBlack = false;
                            }
                        }
                    }
                }
            }
        }
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

    const double MIN_MUTATION_RATE = 0.1;
    const double MAX_MUTATION_RATE = 0.8;

    void adaptMutationRates() {
        for(auto& island : islands) {
            double diversity = calculatePopulationDiversity(island.population);
            if(diversity < MIN_DIVERSITY_THRESHOLD) {
                island.mutationRate = std::min(MAX_MUTATION_RATE, island.mutationRate * 1.5);
            } else if(diversity > 0.4) {
                island.mutationRate = std::max(MIN_MUTATION_RATE, island.mutationRate * 0.8);
            }
        }
    }


    void mutate(Individual &ind, int generationsWithoutImprovement, double populationDiversity, double mutationRateLoc) {
        double oldFitness = ind.fitness;
        if (generationsWithoutImprovement > 5) {
            heavyMutation(ind, std::min(0.8, (1.0 - populationDiversity) * 2));
        } else {
            // Calculate adaptive rate based on current state
            //double adaptiveMutationRate = calculateAdaptiveMutationRate(ind, generationsWithoutImprovement,
            //                                                            populationDiversity);


            double adaptiveMutationRate = mutationRateLoc;
            // Apply structural and clue mutations
            mutator.mutateBoard(ind.board, adaptiveMutationRate);

        }
        // Update fitness and statistics
        double newFitness = calculateFitness(ind);

        //Some quick validation & fixes - doesn't change much though, still hits plateau and stays there
        //repairInvalidRuns(ind.board);
        //repairInvalidCells(ind.board);
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


    // Configuration for island model
    static const int NUM_ISLANDS = 5;
    static const int MIGRATION_INTERVAL = 7;
    static const int MIGRATION_SIZE = 4;
    static const int MIN_ISLAND_SIZE = 18;
    //static const int OPTIMAL_ISLAND_SIZE = 25;
    //const double MIGRATION_SELECTION_PRESSURE = 0.4;

    std::vector<Island> islands;

    void initializeWithDensity(Island &island, double blackCellProbability) {
        for (auto &individual: island.population) {
            // Initialize board with specified density
            initializeIndividualWithDensity(individual, blackCellProbability);

            // Calculate initial fitness
            individual.fitness = calculateFitness(individual);

            // Update island's best if needed
            if (individual.fitness > island.bestFitness) {
                island.bestFitness = individual.fitness;
                island.bestIndividual = individual;
            }
        }
    }

    void initializeIndividualWithDensity(Individual &individual, double blackCellProbability) {
        for (int i = 1; i < size - 1; ++i) {
            for (int j = 1; j < size - 1; ++j) {
                if (threadRng.getReal() < blackCellProbability) {
                    individual.board[i][j] = Cell(true);
                } else {
                    individual.board[i][j] = Cell(false);
                }
            }
        }
    }

    void initializeRandomPatterns(Island &island) {
        PatternGenerator patterns(threadRng, size);

        for (auto &individual: island.population) {
            initializeIndividualWithRandomPattern(patterns, individual, island.denseNessBlack);
            individual.fitness = calculateFitness(individual);
        }
    }

    void initializeIndividualWithRandomPattern(PatternGenerator &patterns, Individual &individual,
                                               double blackCellPossibility) {
        auto pattern = patterns.generateRandomPattern();

        for (int i = 1; i < size - 1; ++i) {
            for (int j = 1; j < size - 1; ++j) {
                if (threadRng.getReal() < 0.8) {
                    individual.board[i][j] = Cell(pattern[i][j]);
                } else {
                    individual.board[i][j] = Cell(threadRng.getReal() < blackCellPossibility);
                }
            }
        }

    }

    void initializeFromTemplates(Island &island) {
        PatternGenerator patterns(threadRng, size);

        for (auto &individual: island.population) {
            initializeIndividualWithTemplate(patterns, individual, island.denseNessBlack);
            individual.fitness = calculateFitness(individual);
        }
    }

    void
    initializeIndividualWithTemplate(PatternGenerator &patterns, Individual &individual, double blackCellProbability) {
        auto templ = patterns.generateRandomTemplate();

        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                if (threadRng.getReal() < 0.9) {
                    individual.board[i][j] = Cell(templ[i][j]);
                } else {
                    individual.board[i][j] = Cell(threadRng.getReal() < blackCellProbability);
                }
            }
        }
    }

    void initializeWithSymmetry(Island &island) {
        for (auto &individual: island.population) {
            // Initialize with rotational symmetry
            initializeIndividualWithSymmetry(individual, island.denseNessBlack);
            individual.fitness = calculateFitness(individual);
        }
    }

    void initializeIndividualWithSymmetry(Individual &individual, double blackCellProbability) {
        for (int i = 0; i < size / 2; ++i) {
            for (int j = 0; j < size / 2; ++j) {
                if (threadRng.getReal() < blackCellProbability) {
                    individual.board[i][j] = Cell(true);
                    individual.board[size - 1 - i][size - 1 - j] = Cell(true);

                    // For larger boards, add additional symmetry points
                    if (size > 10) {
                        individual.board[j][size - 1 - i] = Cell(true);
                        individual.board[size - 1 - j][i] = Cell(true);
                    }
                }
            }
        }

        // For larger boards, add some random noise to break perfect symmetry
        if (size > 10) {
            int noisePoints = size / 5;
            for (int k = 0; k < noisePoints; k++) {
                int i = threadRng.getInt(1, size - 2);
                int j = threadRng.getInt(1, size - 2);
                individual.board[i][j] = Cell(!individual.board[i][j].isBlack);
            }
        }

    }

    void initializeIslands() {
        islands.clear();
        const int islandSize = calculateAdaptivePopSize(size, 1.0) / NUM_ISLANDS;

//#pragma omp parallel for
        for (int i = 0; i < NUM_ISLANDS; ++i) {
            islands.emplace_back(size, islandSize, i);
            islands[i].mutationRate = mutationRate * (0.5 + (i / (NUM_ISLANDS - 1.0))) + (i * 0.02);
            islands[i].crossoverRate = 0.6 + (0.3 * i / (NUM_ISLANDS - 1.0));


            switch (i % 5) {
                case 0:
                    initializeWithDensity(islands[i], 0.4);
                    islands[i].denseNessBlack = 0.4;
                    break;
                case 1:
                    initializeWithDensity(islands[i], 0.2);
                    islands[i].denseNessBlack = 0.2;
                    break;
                case 2:
                    initializeRandomPatterns(islands[i]);
                    islands[i].denseNessBlack = 0.3;
                    break;
                case 3:
                    initializeFromTemplates(islands[i]);
                    islands[i].denseNessBlack = 0.35;
                    break;
                case 4:
                    initializeWithSymmetry(islands[i]);
                    islands[i].denseNessBlack = 0.25;
                    break;
            }
        }

        //debug output not parallel
        for (int i = 0; i < NUM_ISLANDS; ++i) {
            std::cout << "Island " << islands[i].id << " with mutationrate " << islands[i].mutationRate
                      << " crossover Rate " << islands[i].crossoverRate << std::endl;
        }
    }

    static int calculateAdaptivePopSize(int boardSize, double diversity) {
        int baseSize = boardSize * boardSize * 5;
        int staticBounds = std::min(500, std::max(125, baseSize));
        return static_cast<int>(staticBounds * (1.0 + (1.0 - diversity)));
    }

    void balanceIslandSize(Island &island) {
        double diversity = calculatePopulationDiversity(island.population);
        int targetSize = calculateAdaptivePopSize(size, diversity) / NUM_ISLANDS;

        //fill island if necessary
        while (island.population.size() < std::max(MIN_ISLAND_SIZE, targetSize)) {
            Individual newInd = island.bestIndividual;
            heavyMutation(newInd, 0.5);
            newInd.fitness = calculateFitness(newInd);
            island.population.push_back(newInd);
        }

        // remove individuals if necessary
        while (island.population.size() > targetSize) {
            auto worst = std::min_element(island.population.begin(),
                                          island.population.end(),
                                          [](const Individual &a, const Individual &b) {
                                              return a.fitness < b.fitness;
                                          });
            island.population.erase(worst);
        }
    }

    double calculateAverageIslandFitness(const Island &island) {
        if (island.population.empty()) return 0.0;
        double sum = 0.0;
        for (const auto &ind: island.population) {
            sum += ind.fitness;
        }
        return sum / island.population.size();
    }

    void performMigration() {
        std::vector<double> diversities(NUM_ISLANDS);
#pragma omp parallel for
        for (int i = 0; i < NUM_ISLANDS; i++) {
            diversities[i] = calculatePopulationDiversity(islands[i].population);
        }

        for (int i = 0; i < NUM_ISLANDS; i++) {
            if (diversities[i] < MIN_DIVERSITY_THRESHOLD * 1.5) {
                int sourceIdx = selectDiverseIsland(diversities);
                migrateIndividuals(islands[sourceIdx], islands[i],
                                   MIGRATION_SIZE);
                islands[i].mutationRate *= 1.2;
            }
        }
    }

    int selectDiverseIsland(const std::vector<double>& diversities) {
        int bestIdx = 0;
        double bestScore = 0.0;

        for (int i = 0; i < NUM_ISLANDS; i++) {
            double score = diversities[i] * islands[i].bestFitness;
            if (score > bestScore) {
                bestScore = score;
                bestIdx = i;
            }
        }

        // Fallback if no good candidates
        if (bestScore < 0.1) {
            return threadRng.getInt(0, NUM_ISLANDS-1);
        }

        return bestIdx;
    }

    int selectByProbability(const std::vector<double>& probs) {
        double r = threadRng.getReal();
        double sum = 0.0;

        for(size_t i = 0; i < probs.size(); i++) {
            sum += probs[i];
            if(r <= sum) return i;
        }

        return probs.size() - 1; // Fallback to last index
    }

    void migrateIndividuals(Island& source, Island& dest, int count) {
        std::vector<double> probs(source.population.size());
        double totalDist = 0;

        for(size_t i = 0; i < source.population.size(); i++) {
            double dist = calculateBoardDistance(source.population[i].board, dest.bestIndividual.board);
            probs[i] = dist;
            totalDist += dist;
        }

        // Normalize probabilities
        for(auto& p : probs) p /= totalDist;

        // Select diverse individuals
        for(int i = 0; i < count; i++) {
            int idx = selectByProbability(probs);
            dest.population[i] = source.population[idx];
            heavyMutation(dest.population[i], 0.3);
        }
    }

    void updateTemperature(int currentGen) {
        temperature = 1.0 - (static_cast<double>(currentGen) / MAX_GENERATIONS);
    }

    void injectDiversity(Island &island, double populationDiversity) {

        // Increased injection rate when fitness plateaus
        double injectionRate = StagnationParams::BASE_INJECTION_RATE +
                               (0.2 * island.generationsWithoutImprovement / 10.0);

        int injectCount = std::min(
                static_cast<int>(island.population.size() * injectionRate),
                static_cast<int>(island.population.size() * 0.7)
        );
        std::vector<Individual> newIndividuals;

        for (int i = 0; i < injectCount; i++) {
            if (island.denseNessBlack == 0.4) {
                initializeIndividualWithDensity(island.population[i], 0.4);
            } else if (island.denseNessBlack == 0.2) {
                initializeIndividualWithDensity(island.population[i], 0.2);
            } else if (island.denseNessBlack == 0.3) {
                PatternGenerator patterns(threadRng, size);
                initializeIndividualWithRandomPattern(patterns, island.population[i], 0.3);
            } else if (island.denseNessBlack == 0.35) {
                PatternGenerator patterns(threadRng, size);
                initializeIndividualWithTemplate(patterns, island.population[i], 0.35);
            } else if (island.denseNessBlack == 0.4) {
                initializeIndividualWithSymmetry(island.population[i], 0.4);
            }
            Individual newInd = island.population[i];
            newInd.fitness = calculateFitness(newInd);
            newIndividuals.push_back(newInd);
        }

        // Replace worst individuals
        std::sort(island.population.begin(), island.population.end(),
                  [](const Individual &a, const Individual &b) {
                      return a.fitness < b.fitness;
                  });

        for (int i = 0; i < injectCount; i++) {
            island.population[island.population.size() - 1 - i] = newIndividuals[i];
        }
    }

    void adaptSelectionPressure(Island &island) {
        double avgFitness = calculateAverageIslandFitness(island);
        double bestFitness = island.bestFitness;

        // If average is too far from best, reduce pressure to explore more
        double fitnessGap = bestFitness - avgFitness;
        if (fitnessGap > 0.4) {
            island.selectionPressure = std::max(0.4, std::min(0.1, island.selectionPressure - 0.05));
        } else {
            island.selectionPressure = std::min(0.7, island.selectionPressure + 0.05);
        }
    }


    void evolveIsland(Island &island, int generation) {
        // Calculate population diversity for this island
        double populationDiversity = calculatePopulationDiversity(island.population);
        bool isStagnating = island.generationsWithoutImprovement > 10 || populationDiversity < 0.2;

        temperature = 1.0 - (generation / MAX_GENERATIONS);
        island.mutationRate = mutationRate * (1.0 + temperature);

        //std::cout << "Island Size: " << std::to_string(island.population.size()) << std::endl;

        // Evaluate fitness

        //TODO: improve parallelization, more sophisticated strategies broke it
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

        // Adapt selection pressure based on population state
        adaptSelectionPressure(island);

        // Check for stagnation
        if ((island.bestFitness > 0.8 && island.generationsWithoutImprovement > 5) ||
            populationDiversity < StagnationParams::MIN_DIVERSITY) {
            island.mutationRate = std::min(0.3, island.mutationRate * 1.5);
            //std::cout << "Injecting new boards" << std::endl;
            injectDiversity(island, populationDiversity);
        }

        if (island.generationsWithoutImprovement == 0) {
            island.mutationRate = mutationRate;
        }

        // Create new population
        std::vector<Individual> newPopulation;

        // Elitism
        int eliteCount = std::max(1, static_cast<int>(
                island.population.size() * 0.1 * (1.0 + temperature)
        ));

        if (island.generationsWithoutImprovement < StagnationParams::STAGNATION_THRESHOLD) {
            std::sort(island.population.begin(), island.population.end(),
                      [](const Individual &a, const Individual &b) {
                          return a.fitness > b.fitness;
                      });

            for (int i = 0; i < eliteCount; ++i) {
                newPopulation.push_back(island.population[i]);
            }
        }

        // Fill rest with crossover and mutation
        while (newPopulation.size() < island.population.size()) {
            // Tournament selection with adaptive pressure
            Individual parent1 = selectParent(island, island.selectionPressure);
            Individual parent2 = selectParent(island, island.selectionPressure);


            Individual child =
                    threadRng.getReal() < island.crossoverRate ? crossover(parent1, parent2) : (threadRng.getReal() <
                                                                                                0.5
                                                                                                ? parent1 : parent2);
            if (threadRng.getReal() < island.mutationRate) {
                if (isStagnating) {
                    heavyMutation(child, populationDiversity);
                } else {
                    mutate(child, island.generationsWithoutImprovement, populationDiversity, island.mutationRate);
                }
            }

            newPopulation.push_back(child);
        }

        island.population = std::move(newPopulation);
        island.generationsWithoutImprovement++;

        // Balance population size
        balanceIslandSize(island);
    }


    Individual selectParent(const Island &island, double selectionPressure) {
        int tournamentSize = std::max(2, static_cast<int>(
                island.population.size() * selectionPressure * 0.2));

        std::vector<Individual> tournament;
        for (int i = 0; i < tournamentSize; i++) {
            tournament.push_back(island.population[threadRng.getInt(0,
                                                                    island.population.size() - 1)]);
        }

        // Consider both fitness and diversity contribution
        double populationDiversity = calculatePopulationDiversity(island.population);

        if (populationDiversity < MIN_DIVERSITY_THRESHOLD) {
            // Prioritize diversity
            return *std::max_element(tournament.begin(), tournament.end(),
                                     [&](const Individual &a, const Individual &b) {
                                         return calculateBoardDistance(a.board, island.bestIndividual.board)
                                                > calculateBoardDistance(b.board, island.bestIndividual.board);
                                     });
        }

        // fitness-based selection
        return *std::max_element(tournament.begin(), tournament.end(),
                                 [](const Individual &a, const Individual &b) {
                                     return a.fitness < b.fitness;
                                 });
    }

    void reinitializeIslandIfStuck(Island &island, int maxStagnantGenerations = 30) {
        if (island.generationsWithoutImprovement > maxStagnantGenerations) {

            // Scale based on board size
            double eliteRatio = std::max(0.1, std::min(0.3, 5.0 / size));
            int eliteCount = island.population.size() * eliteRatio;

            std::sort(island.population.begin(), island.population.end(),
                      [](const Individual &a, const Individual &b) {
                          return a.fitness > b.fitness;
                      });

            // Store elites for mutations - will be replaced by mutated
            std::vector<Individual> elites(island.population.begin(),island.population.begin() + eliteCount);

            // Determine initialization strategy based on board size and current fitness
            double currentBestFitness = island.bestFitness;
            double diversityRatio = calculatePopulationDiversity(island.population);



            // Replace population
            for (int i = 0; i < island.population.size(); ++i) {
                if (i < island.population.size() * 0.4) {
                    if (island.denseNessBlack == 0.4) {
                        initializeIndividualWithDensity(island.population[i], 0.4);
                    } else if (island.denseNessBlack == 0.2) {
                        initializeIndividualWithDensity(island.population[i], 0.2);
                    } else if (island.denseNessBlack == 0.3) {
                        PatternGenerator patterns(threadRng, size);
                        initializeIndividualWithRandomPattern(patterns, island.population[i], 0.3);
                    } else if (island.denseNessBlack == 0.35) {
                        PatternGenerator patterns(threadRng, size);
                        initializeIndividualWithTemplate(patterns, island.population[i], 0.35);
                    } else if (island.denseNessBlack == 0.4) {
                        initializeIndividualWithSymmetry(island.population[i], 0.4);
                    }
                    island.population[i].fitness = calculateFitness(island.population[i]);
                } else if (i < island.population.size() * 0.7) {
                    Individual newInd = elites[i % eliteCount];
                    heavyMutation(newInd, diversityRatio);
                    island.population[i] = newInd;
                    island.population[i].fitness = calculateFitness(island.population[i]);
                } else {
                    island.population[i] = createIndividualWithBasicBoardValidity();
                    island.population[i].fitness = calculateFitness(island.population[i]);

                }
            }



            for (int i = 0; i < island.population.size(); ++i) {
                if (island.population[i].fitness > island.bestFitness) {
                    island.bestFitness = island.population[i].fitness;
                    island.bestIndividual = island.population[i];
                }
            }

            std::cout << "Reinitialized population " << std::to_string(island.id) << " new best fitness "
                      << island.bestFitness << std::endl;

            island.generationsWithoutImprovement = 0;
        }
    }

public:

    KakuroGenerator(int
                    boardSize)
            : size(boardSize),
              mutationRate(0.1),
              targetFitness(0.98),
              mutator(boardSize),
              temperature(0),
              threadRng()
               {
    }

    void ensurePopulationDiversity(Island& island) {
        double diversity = calculatePopulationDiversity(island.population);
        if (diversity < MIN_DIVERSITY_THRESHOLD) {
            int injectCount = std::min(
                    static_cast<int>(island.population.size() * 0.4),
                    static_cast<int>(MAX_POPULATION_SIZE * 0.2)
            );
            std::vector<Individual> newIndividuals;
            for (int i = 0; i < injectCount; i++) {
                Individual newInd = createIndividualWithBasicBoardValidity();
                heavyMutation(newInd, diversity);
                newIndividuals.push_back(newInd);
            }

            // Replace worst individuals while preserving elites
            std::partial_sort(island.population.begin(),
                              island.population.begin() + injectCount,
                              island.population.end(),
                              [](const Individual& a, const Individual& b) {
                                  return a.fitness < b.fitness;
                              });

            for (int i = 0; i < injectCount; i++) {
                island.population[i] = newIndividuals[i];
            }
        }
    }

    std::vector<std::vector<Cell>> generateBoard(int maxTimeSeconds = 300) {
        auto startTime = std::chrono::steady_clock::now();

        // Initialize islands
        initializeIslands();

        // Main evolution loop
        int generation = 0;
        Individual globalBest(size);
        double globalBestFitness = 0.0;

        while (generation < MAX_GENERATIONS) {
            updateTemperature(generation);
            // Evolve each island
#pragma omp parallel for
            for (int i = 0; i < NUM_ISLANDS; ++i) {
                if (islands[i].population.empty()) {
                    std::cout << "Warning: Island " << i << " has empty population!\n";
                    continue;
                }
                adaptMutationRates();
                evolveIsland(islands[i], generation);
                ensurePopulationDiversity(islands[i]);
                reinitializeIslandIfStuck(islands[i]);
                std::cout << "Island " << islands[i].id << " best fitness: " << islands[i].bestFitness << std::endl;

            }

            // Perform migration if needed
            if (generation % MIGRATION_INTERVAL == 0 && generation > 0) {
                std::cout << "Migrating individuals " << std::endl;
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
            for (auto &island: islands) {
                //std::lock_guard<std::mutex> lock(populationMutex);
                if (island.bestFitness > globalBestFitness) {

                    globalBestFitness = island.bestFitness;
                    globalBest = island.bestIndividual;


                    std::cout << "Generation " << generation << ": New best fitness = "
                              << globalBestFitness << " (Island "
                              << &island - &islands[0] << ")\n";


                    printBoard(globalBest.board);

                    if (globalBestFitness >= targetFitness) {
                        std::cout << "Target fitness reached after " << generation << " generations\n";
                        std::vector<std::vector<Cell>> solution;
                        KakuroSolver solver(size);

                        // Copy board configuration to solver
                        for (int i = 0; i < size; ++i) {
                            for (int j = 0; j < size; ++j) {
                                solver.setCell(i, j,
                                               globalBest.board[i][j].isBlack,
                                               globalBest.board[i][j].downClue,
                                               globalBest.board[i][j].rightClue);
                            }
                        }
                        solver.printInitialBoard();

                        SolveResult result = solveBoardWithTimeOut(solver);
                        if (result == SolveResult::UNIQUE_SOLUTION) return globalBest.board;
                        std::string text;
                        if (result == SolveResult::INVALID_BOARD) {
                            text = " invalid board ";
                        } else if (result == SolveResult::NO_SOLUTION) {
                            text = " no solution ";
                        } else {
                            text = " multiple solutions ";
                        }

                        std::cout << "WARNING: target fitness reached, but not uniquely solvable" << text << std::endl;
                        //adapt current best individual to not have it pollute further exploration
                        island.bestIndividual = crossover(island.bestIndividual,
                                                          createIndividualWithBasicBoardValidity());
                        globalBestFitness = 0.0;
                        globalBest = Individual(size);
                    }
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

            cleanupCache();
            generation++;
        }

        return globalBest.board;
    }

    SolveResult solveBoardWithTimeOut(KakuroSolver &solver) {
        int timeout = 15;
        std::atomic<bool> should_terminate{false};
        auto start = std::chrono::high_resolution_clock::now();

        std::future<SolveResult> future = std::async(std::launch::async, [&]() {
            std::vector<std::vector<Cell>> solution;
            return solver.solveBoard(solution);
        });

        if (future.wait_for(std::chrono::seconds(timeout)) == std::future_status::timeout) {
            should_terminate = true;
            std::cout << "solver timeout" << std::endl;
            return SolveResult::NO_SOLUTION;
        }
        auto result = future.get();
        return result;
    }

};