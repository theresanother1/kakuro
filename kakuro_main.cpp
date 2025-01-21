#include <iostream>
#include <vector>
#include <string>
#include "KakuroSolver.cpp"
#include "KakuroGenerator.cpp"

/*
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
*/


void solveCreatedKakuroCheck(std::string &filename, std::string &solutionFilename) {
    auto solver = KakuroSolver::readFromFile(filename);

    std::cout << "\nInitial board setup:\n";
    solver.printInitialBoard();

    std::vector<std::vector<Cell>> solution;
    SolveResult hasUniqueSolution = solver.solveBoard(solution);

    if (hasUniqueSolution == SolveResult::UNIQUE_SOLUTION) {
        std::cout << "\nUnique solution found!\n";
        //solver.writeToFile(solutionFilename, solution);
    } else if (hasUniqueSolution == SolveResult::MULTIPLE_SOLUTIONS) {
        std::cout << "\nMultiple solutions exist.\n";
    } else {
        std::cout << "\nNo solution exists.\n";
    }

}


int main() {
    std::string file = "old_files/initial_test_6x6.ka";
    std::string second = "old_files/board_5x5.kakuro";

    solveCreatedKakuroCheck(file, file);
    solveCreatedKakuroCheck(second, file);

    int run = 2;
    SolveResult result = SolveResult::NO_SOLUTION;

    for (int size = 5; size <= 20; ++size) {

        KakuroBoard boardTest(size, size);
        boardTest.generateBoard();
        std::cout << "Starting generating board size " << std::to_string(size) << std::endl;
        do {
            KakuroGenerator generator(size);
            auto board = generator.generateBoard();

            std::cout  << "Finished generating, try result" << std::endl;

            KakuroSolver solver(size);
            for (int i = 0; i < size; i++) {
                for (int j = 0; j < size; j++) {
                    solver.setCell(i, j, board[i][j].isBlack,
                                   board[i][j].downClue, board[i][j].rightClue);
                }
            }
            std::vector <std::vector<Cell>> solution;
            result = solver.solveBoard(solution);
            if (result == SolveResult::UNIQUE_SOLUTION) {
                std::string fileName = std::to_string(run)+"_board_" + std::to_string(size) +"x"+std::to_string(size)+".kakuro";
                std::string fileNameSol = std::to_string(run)+"_board_" + std::to_string(size) +"x"+std::to_string(size)+"_solution.kakuro";
                solver.writeToFile(fileName, board);
                solver.writeToFile(fileNameSol, solution);
            }
        } while (result != SolveResult::UNIQUE_SOLUTION);
        std::cout << "Finished generating board size " << std::to_string(size) << std::endl;
    }

    return 0;
}