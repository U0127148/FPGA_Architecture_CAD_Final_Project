#include "solver.h"

int main(int argc, char *argv[]) {
    srand(0);
    Solver* solver = new Solver();
    solver->solve(argv);
    return 0;
}