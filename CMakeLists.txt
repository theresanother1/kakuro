cmake_minimum_required(VERSION 3.28)
project(part2)

set(CMAKE_CXX_STANDARD 17)

include_directories(.)

add_executable(part2
        kakuro_main.cpp
        KakuroSolver.h
        KakuroSolver.cpp
        KakuroGenerator.cpp
)
