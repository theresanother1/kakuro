Makefile included --> should work on Windows and Unix based. 

kakuro_main.cpp includes main function: 
* executes solving two boards, first solvable, second isn't to show solver capabilities
* then executes evolutionary loops for each size with output to console 

Exercise description:
Please create kakuro boards (size 5x5  up to 20x20) with an EA that can be solved in a unique way (so that there is only one solution) in C++.

Can be executed two ways: 
- use the executable only --> executes the evolutionary algorithm in a loop for all dimensions
- provide a filename --> tries to solve the given puzzle 