#include <iostream>
#include <vector>
#include <string>
#include "KakuroSolver.cpp"
#include "KakuroGenerator.cpp"

void solveCreatedKakuroCheck(std::string &filename, std::string &solutionFilename) {
    KakuroSolver solver = KakuroSolver::readFromFile(filename);

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


int main(int argc, char *argv[]) {
    try {
        if (argc == 1) {
            std::string file = "old_files/initial_test_6x6.ka";
            std::string second = "old_files/board_5x5.kakuro";

            solveCreatedKakuroCheck(file, file);
            solveCreatedKakuroCheck(second, file);


            int run = 2;
            SolveResult result = SolveResult::NO_SOLUTION;

            for (int size = 5; size <= 20; ++size) {
                KakuroBoard boardTest(size);
                boardTest.generateBoard();
                std::cout << "Starting generating board size " << std::to_string(size) << std::endl;
                do {
                    KakuroGenerator generator(size);
                    auto board = generator.generateBoard();

                    std::cout << "Finished generating, try result" << std::endl;

                    KakuroSolver solver(size);
                    std::vector<std::vector<Cell>> solution;
                    solver.initializeBoard(board);
                    result = solver.solveBoard(solution);
                    if (result == SolveResult::UNIQUE_SOLUTION) {
                        std::string fileName =
                                std::to_string(run) + "_board_" + std::to_string(size) + "x" + std::to_string(size) +
                                ".kakuro";
                        std::string fileNameSol =
                                std::to_string(run) + "_board_" + std::to_string(size) + "x" + std::to_string(size) +
                                "_solution.kakuro";
                        solver.writeToFile(fileName, board);
                        solver.writeToFile(fileNameSol, solution);
                    }
                    std::cout << "Finished generating board size " << std::to_string(size) << std::endl;
                } while (result != SolveResult::UNIQUE_SOLUTION);
                std::cout << "Next Board " << std::endl;
            }
        } else {
            std::string check = argv[1];
            //std::cin >> check;
            std::string fileE = " ";
            solveCreatedKakuroCheck(check, fileE);
        }
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
    }


    return 0;
}