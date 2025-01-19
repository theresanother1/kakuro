#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include "KakuroSolver.cpp"
#include "KakuroGenerator.cpp"

void solveExampleKakuro() {
    KakuroSolver solver(5);

    // First set all cells as black (done in constructor)

    // Set white cells
    solver.setCell(1, 2, false);  // Row 1
    solver.setCell(1, 3, false);
    solver.setCell(1, 4, false);

    solver.setCell(2, 2, false);  // Row 2
    solver.setCell(2, 3, false);
    solver.setCell(2, 4, false);

    solver.setCell(3, 1, false);  // Row 3
    solver.setCell(3, 2, false);
    solver.setCell(3, 3, false);

    solver.setCell(4, 1, false); // Row 4
    solver.setCell(4, 2, false);
    solver.setCell(4, 3, false);

    // Set clues - verified solution exists
    solver.setCell(0, 2, true, 10, 0);    // Down clue 10
    solver.setCell(0, 3, true, 30, 0);    // Down clue 30
    solver.setCell(0, 4, true, 4, 0);    // Down clue 4

    solver.setCell(1, 1, true, 0, 10);    // Right clue 10
    solver.setCell(2, 1, true, 17, 15);    // Right clue 15 / Down 17
    solver.setCell(3, 0, true, 0, 19);    // Right clue 19
    solver.setCell(4, 0, true, 0, 17);     // Right clue 17

    std::cout << "\nInitial board setup:\n";
    // solver.printInitialBoard();

    std::vector<std::vector<Cell>> solution;
    bool hasUniqueSolution = solver.solveBoard(solution);

    if (hasUniqueSolution) {
        std::cout << "\nUnique solution found!\n";
    } else {
        std::cout << "\nNo unique solution exists.\n";
    }

}

void solveCreatedKakuroCheck(std::string &filename, std::string &solutionFilename) {
    auto solver = KakuroSolver::readFromFile(filename);

    std::cout << "\nInitial board setup:\n";
    solver.printInitialBoard();

    std::vector<std::vector<Cell>> solution;
    bool hasUniqueSolution = solver.solveBoard(solution);

    if (hasUniqueSolution) {
        std::cout << "\nUnique solution found!\n";
        solver.writeToFile(solutionFilename, solution);
    } else {
        std::cout << "\nNo unique solution exists.\n";
    }

}


void generateBoardMap() {
    std::cout << "Starting Kakuro board map generation...\n\n";
    int run = 5;

    // For each size from 5x5 to 20x20
    for (int size = 6; size <= 20; ++size) {
        std::string filename = std::to_string(run) + "_board_" + std::to_string(size) + "x" +
                               std::to_string(size) + ".kakuro";

        std::string filename_solution = std::to_string(run) + "_board_" + std::to_string(size) + "x" +
                                        std::to_string(size) + "solution.kakuro";

        std::cout << "Generating " << size << "x" << size << " board...\n";

        // Adjust population size based on board size
        int populationSize = 50 + (size - 5) * 10;  // Larger populations for bigger boards
        KakuroGenerator generator(size, populationSize);

        auto start = std::chrono::high_resolution_clock::now();
        // Generate board with increased time limit for larger boards
        //auto timeLimit = std::chrono::seconds(300 + (size - 5) * 60);  // 5 mins + 1 min per size above 5
        auto board = generator.generateBoard();

        std::cout << "Created board with " << size << "x" << size << std::endl;
        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Time taken: " << duration.count() << " ms\n";


        generator.writeToFile(filename, board);

        // Create solver to verify the board
        KakuroSolver solver(size);
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                solver.setCell(i, j, board[i][j].isBlack,
                               board[i][j].downClue, board[i][j].rightClue);
            }
        }

        // Verify and save board
        std::vector<std::vector<Cell>> solution;


        start = std::chrono::high_resolution_clock::now();
        bool isValid = solver.solveBoard(solution);

        end = std::chrono::high_resolution_clock::now();

        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Time taken: " << duration.count() << " ms\n";

        if (isValid) {
            solver.writeToFile(filename_solution, solution);
            std::cout << "Successfully generated and saved " << filename_solution << "\n\n";
        } else {
            std::cerr << "Failed to generate valid " << size << "x" << size
                      << " board. \n";

        }
    }

    std::cout << "Board map generation completed!\n";
}


int main() {
    std::cout << "Check - 1 or run algo - 0 ";
    int check = 0;
    //std::cin >> check;

    if (check == 1) {
        auto start = std::chrono::high_resolution_clock::now();
        solveExampleKakuro();
        auto end = std::chrono::high_resolution_clock::now();


        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Time taken: " << duration.count() << " ms\n";

        start = std::chrono::high_resolution_clock::now();
        std::string filename = "1_board_5x5.kakuro";
        std::string filenameSolution = "1_board_5x5solution.kakuro";
        std::string filenameWrong = "wrong_board.txt";
        std::string filenameWrongSol = "wrong_board_solution.txt";



        solveCreatedKakuroCheck(filename, filenameSolution);
        solveCreatedKakuroCheck(filenameWrong, filenameWrongSol);
        end = std::chrono::high_resolution_clock::now();

        duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Time taken: " << duration.count() << " ms\n";


    } else {
        generateBoardMap();
    }

    return 0;
}
