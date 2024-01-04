#include <iostream>
#include <string>
#include <vector>
#include <climits>

enum blockType {
    IO,
    CLB,
    RAM,
    DSP
};

struct Block {
    std::string name;
    int type;
    double center_x, center_y;
    int id;
    int sel;
    Block(std::string n, int t, double x, double y) : name(n), type(t), center_x(x), center_y(y), sel(0) {}
};

struct Net {
    std::string name;
    std::vector<Block*> inst_vec;
    Net(std::string name) : name(name) {}
};

struct gene {
    double fitness;
    std::vector<int> resource_permu[4]; // store resource block (vector index means the instance)
};

typedef std::pair<gene, gene> parents;