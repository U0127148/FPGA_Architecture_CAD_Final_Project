RES = solver.cpp main.cpp
EXE = legalizer
all:
	g++ -O3 -std=c++11 $(RES) -o $(EXE)
clean:
	rm $(EXE)