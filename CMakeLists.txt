cmake_minimum_required(VERSION 3.28)
project(part2)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(part2
        kakuro_main.cpp
        KakuroSolver.h
        KakuroSolver.cpp
        KakuroGenerator.h
        KakuroGenerator.cpp
)
