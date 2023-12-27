#include "solver.h"

int main(int argc, char *argv[]) {
    srand(0);
    Solver* solver = new Solver();
    solver->read(argv);
    solver->create_grid();
    solver->init_pop_from_GP();
    solver->output_file(argv[4]);
    solver->genetic_algorithm();
    solver->output_file(argv[4]);
    return 0;
}
