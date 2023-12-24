#include "solver.h"

int main(int argc, char *argv[]) {
    srand(0);
    Solver* solver = new Solver();
    solver->read(argv);
    solver->init_pop();
    solver->genetic_algorithm();
    solver->output_file(argv[4]);
    return 0;
}