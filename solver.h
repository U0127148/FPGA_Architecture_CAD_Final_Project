#include "data.h"
#include <unordered_map>

class Solver {
public:
    Solver() {}
    ~Solver() {}
    void read(char *argv[]);
    void makeWindow();
    void setUpObject();
private:
    std::vector<Block*> resource;
    std::vector<Block*> inst;
    std::unordered_map<std::string, Block*> nameToResource;
    std::unordered_map<std::string, Block*> nameToInst;
    std::vector<Net*> net_vec;

    std::vector<std::vector<Block*>> window[4];

    int mcell_num;
    int **C_matrix;
    int **d_x, **d_y;
};