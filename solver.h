#include "data.h"
#include <unordered_map>
#include <chrono>

class Solver {
public:
    Solver() {}
    ~Solver() {}
    void read(char *argv[]);
    void create_grid();
    int find_block(int, int, int); // type, inst_idx, pop_size
    int find_available_block(int, int, int, int); // type row col pop_size
    void init_pop();
    void fitness(gene&);
    bool incremental_fitness(gene&, int, int, int); // type, instance index, resource index
    void parent_selection(parents&);
    void crossover(parents&, std::vector<gene>&);
    void mutation(parents&, std::vector<gene>&);
    void survivor_selection(std::vector<gene>&);
    void genetic_algorithm(double);
    void output_file(char*);
    void solve(char *argv[]);
private:
    std::vector<Block*> resource[4];
    std::vector<Block*> inst[4];
    std::unordered_map<std::string, Block*> nameToResource;
    std::unordered_map<std::string, Block*> nameToInst;
    std::vector<Net*> net_vec;

    std::vector<std::vector<int>> inst_to_net[4];

    std::vector<gene> pool;
};