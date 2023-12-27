#include "data.h"
#include <unordered_map>

class Solver {
public:
    Solver() {}
    ~Solver() {}
    void read(char *argv[]);
    void create_grid();
    int find_block(int, int, int); // type, inst_idx, pop_size
    int find_available_block(int, int, int, int); // type row col pop_size
    void init_pop_from_GP();
    void fitness(gene&);
    void parent_selection(parents&);
    void crossover(parents&, std::vector<gene>&);
    void mutation(std::vector<gene>&);
    void survivor_selection(std::vector<gene>&);
    void genetic_algorithm();
    void output_file(char*);
private:
    std::vector<Block*> resource[4];
    std::vector<Block*> inst[4];
    std::unordered_map<std::string, Block*> nameToResource;
    std::unordered_map<std::string, Block*> nameToInst;
    std::vector<Net*> net_vec;

    std::vector<std::vector<Block*>> window[4];

    std::vector<std::vector<int>> inst_to_net[4];
    std::vector<gene> pool;
};
