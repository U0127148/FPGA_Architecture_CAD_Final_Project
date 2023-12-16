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
    int type;
    double center_x, center_y;
    int layer, unit;
    Block(int t, double x, double y) : type(t), center_x(x), center_y(y) {}
};

struct Net {
    std::string name;
    std::vector<Block*> inst_vec;
    double min_x, min_y, max_x, max_y;
    double HPWL;
    Net(std::string name) : name(name), min_x(INT_MAX), min_y(INT_MAX), max_x(INT_MIN), max_y(INT_MIN) {}
    void init() {
        for (auto &inst : inst_vec) {
            if (inst->center_x < min_x) min_x = inst->center_x;
            if (inst->center_x > max_x) max_x = inst->center_x;
            if (inst->center_y < min_y) min_y = inst->center_y;
            if (inst->center_y > max_y) max_y = inst->center_y;
        }
        HPWL = (max_x - min_x) + (max_y - min_y);
    }
};

struct Cluster {
    std::vector<Net*> net_vec;
};