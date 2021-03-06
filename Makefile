#Enables C++11
PATH := /usr/um/gcc-4.8.2/bin:$(PATH)
LD_LIBRARY_PATH := /usr/um/gcc-4.8.2/lib64
LD_RUN_PATH := /usr/um/gcc-4.8.2/lib64

#Default Flags
FLAGS = -Wall -Wextra -pedantic -Wvla -std=c++11

# make release - will compile "all" with $(FLAGS) and the O3 flag
release: FLAGS += -O3
release: all

# make debug - will compile "all" with $(FLAGS) and the g flag
debug: FLAGS += -g
debug: all

all:2D_Full_Euler_CFD.o read_parameters.o read_grid.o Exact_Riemann_Solver.o compute_neighbors.o helper_functions.o write_to_file.o
	g++ $(FLAGS) 2D_Full_Euler_CFD.o read_parameters.o read_grid.o Exact_Riemann_Solver.o compute_neighbors.o helper_functions.o write_to_file.o -o 2D_Full_Euler_CFD

# This is an example from one of my projects. Replace the header files with your
# header files and the cpp files with your cpp files 
2D_Full_Euler_CFD.o: 2D_Full_Euler_CFD.cpp file_header.h helper_functions.h
	g++ $(FLAGS) -c 2D_Full_Euler_CFD.cpp
read_parameters.o: read_parameters.cpp file_header.h helper_functions.h
	g++ $(FLAGS) -c read_parameters.cpp
read_grid.o: read_grid.cpp file_header.h helper_functions.h
	g++ $(FLAGS) -c read_grid.cpp
Exact_Riemann_Solver.o: Exact_Riemann_Solver.cpp file_header.h helper_functions.h
	g++ $(FLAGS) -c Exact_Riemann_Solver.cpp
compute_neighbors.o: compute_neighbors.cpp file_header.h helper_functions.h
	g++ $(FLAGS) -c compute_neighbors.cpp
helper_functions.o: helper_functions.cpp file_header.h helper_functions.h
	g++ $(FLAGS) -c helper_functions.cpp
write_to_file.o: write_to_file.cpp file_header.h helper_functions.h
	g++ $(FLAGS) -c write_to_file.cpp

# ------------------ Test cases! ----------------------------
test: test1
	@echo "\ntests successfully finished!\n"

# this assumes you have a folder called test with input called test.txt and your
# executable's name is executable. Also assumes you have a file called
# test1correct.txt in your test folder that has all the correct outputs.
# You can view your output in test/output.txt
test1: all
	./executable < test/spectest.txt > test/output.txt
	diff test/output.txt test/test1correct.txt

# ------------------ Cleans your folder --------------------
clean:
	rm -f *.o letterman
