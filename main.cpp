#include "solver.h"

int main(int argc, char *argv[]) {
    Solver* solver = new Solver();
    solver->read(argv);
    solver->makeWindow();
    return 0;
}