#include "data.h"
#include <unordered_map>

class Solver {
public:
    Solver() {}
    ~Solver() {}
    void read(char *argv[]);
    void makeWindow();
    void setUpObject();
    void getGlobalMinimum();
private:
    std::vector<Block*> resource;
    std::vector<Block*> inst;
    std::unordered_map<std::string, Block*> nameToResource;
    std::unordered_map<std::string, Block*> nameToInst;
    std::vector<Net*> net_vec;

    std::vector<std::vector<Block*>> window[4];

    int mcell_num = 0, fcell_num = 0;
    double **C_matrix;
    double **d_x, **d_y;
};