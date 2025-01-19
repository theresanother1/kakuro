CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O3 -fopenmp
LDFLAGS = -lgomp

SRC =  kakuro_main.cpp
OBJ = $(SRC:.cpp=.o)
EXEC = kakuro
PROBLEMS_DIR = part2

ifeq ($(OS),Windows_NT)
    RM = del /Q
    EXT = .exe
else
    RM = rm -f
    EXT =
endif

.PHONY: all clean

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CXX) $(OBJ) -o $(EXEC)$(EXT) $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJ) $(EXEC)$(EXT)
