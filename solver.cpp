#include "solver.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <random>
#include <climits>
#include <unordered_set>
#include <unordered_map>
const int POP_SIZE = 100;
const int TERMINATION = 100;
const int K = 2;
double min_y[4];
double height[4];
int unit[4];
int size[4];
std::unordered_set<int> st1, st2;
int *mod_table[4];
std::unordered_map<int, int> m; // (child1 value, child1 index)

void Solver::build_mod_table() {
    int mod = 0;
    // CLB
    mod_table[blockType::CLB] = new int[2 * resource[blockType::CLB].size() + 1];
    for (int val = 0; val <= 2 * resource[blockType::CLB].size(); val++) {
        mod_table[blockType::CLB][val] = mod++;
        if (mod == resource[blockType::CLB].size()) mod = 0;
    }
    // for (int val = 0; val <= 2 * resource[blockType::CLB].size(); val++) std::cout << mod_table[blockType::CLB][val] << ' ';
    // std::cout << std::endl;

    // RAM
    mod = 0;
    mod_table[blockType::RAM] = new int[2 * resource[blockType::RAM].size() + 1];
    for (int val = 0; val <= 2 * resource[blockType::RAM].size(); val++) {
        mod_table[blockType::RAM][val] = mod++;
        if (mod == resource[blockType::RAM].size()) mod = 0;
    }
    // for (int val = 0; val <= 2 * resource[blockType::RAM].size(); val++) std::cout << mod_table[blockType::RAM][val] << ' ';
    // std::cout << std::endl;

    // DSP
    mod = 0;
    mod_table[blockType::DSP] = new int[2 * resource[blockType::DSP].size() + 1];
    for (int val = 0; val <= 2 * resource[blockType::DSP].size(); val++) {
        mod_table[blockType::DSP][val] = mod++;
        if (mod == resource[blockType::DSP].size()) mod = 0;
    }
    // for (int val = 0; val <= 2 * resource[blockType::DSP].size(); val++) std::cout << mod_table[blockType::DSP][val] << ' ';
    // std::cout << std::endl;
}

void Solver::read(char *argv[]) {
    std::string line;

    // first input file (resource)
    std::ifstream fin;
    fin.open(argv[1]);
    while(getline(fin, line)) {
        // std::cout << line << std::endl;
        int n = 0;
        std::string info, name;
        int block_type;
        double center_x, center_y;
        std::stringstream ss(line);
        while(getline(ss, info, ' ')) {
            if (n == 0) {
                // std::cout << "name: " << info << ' ';
                name = info;
            } else if (n == 1) {
                // std::cout << "type: " << info << ' ';
                if (info == "CLB") block_type = blockType::CLB;
                else if (info == "RAM") block_type = blockType::RAM;
                else if (info == "DSP") block_type = blockType::DSP;
                else std::cout << "[ERROR] read FIRST input file\n"; 
            } else if (n == 2) {
                // std::cout << "center x: " << info << ' ';
                center_x = std::stod(info);
            } else {
                // std::cout << "center y: " << info << std::endl;
                center_y = std::stod(info);
                // std::cout << "name: " << name << " type: " << block_type << " center x: " << center_x << " center_y: " << center_y << std::endl;
                Block* resource_block = new Block(name, block_type, center_x, center_y);
                resource[block_type].emplace_back(resource_block);
                nameToResource[name] = resource_block;
            }
            n++;
        }
    }
    fin.close();

    int idx[4] = {0, 0, 0, 0};
    // second input file (inst)
    fin.open(argv[2]);
    while(getline(fin, line)) {
        // std::cout << line << std::endl;
        int n = 0;
        std::string info, name;
        int block_type;
        double center_x, center_y;
        int id;
        std::stringstream ss(line);
        while(getline(ss, info, ' ')) {
            if (n == 0) {
                // std::cout << "name: " << info << ' ';
                name = info;
            } else if (n == 1) {
                // std::cout << "type: " << info << ' ';
                if (info == "CLB") {
                    block_type = blockType::CLB;
                    id = idx[blockType::CLB]++;
                } else if (info == "RAM") {
                    block_type = blockType::RAM;
                    id = idx[blockType::RAM]++;
                } else if (info == "DSP") {
                    block_type = blockType::DSP;
                    id = idx[blockType::DSP]++;
                } else if (info == "IO") {
                    block_type = blockType::IO;
                    id = idx[blockType::IO]++;
                } else {
                    std::cout << "[ERROR] read SECOND input file\n"; 
                }
            } else if (n == 2) {
                // std::cout << "center x: " << info << ' ';
                center_x = std::stod(info);
            } else {
                // std::cout << "center y: " << info << std::endl;
                center_y = std::stod(info);
                // std::cout << "name: " << name << " type: " << block_type << " center x: " << center_x << " center_y: " << center_y << std::endl;
                Block* inst_block = new Block(name, block_type, center_x, center_y);
                inst[block_type].emplace_back(inst_block);
                nameToInst[name] = inst_block;
                inst_block->id = id;
            }
            n++;
        }
    }
    fin.close();

    // for (auto& cell : inst[blockType::CLB]) std::cout << cell->id << ' ';
    // std::cout << std::endl;
    // for (auto& cell : inst[blockType::RAM]) std::cout << cell->id << ' ';
    // std::cout << std::endl;
    // for (auto& cell : inst[blockType::DSP]) std::cout << cell->id << ' ';
    // std::cout << std::endl;

    // third input file (net info)
    fin.open(argv[3]);
    Net* net;
    while(getline(fin, line)) {
        // std::cout << line << std::endl;
        int n = 0;
        std::string info, name;
        std::stringstream ss(line);
        while(getline(ss, info, ' ')) {
            if (n == 0) {
                net = new Net(info);
            } else {
                net->inst_vec.emplace_back(nameToInst[info]);
            }
            n++;
        }
        net_vec.emplace_back(net);
    }
    fin.close();

    // CLB
    inst_to_net[blockType::CLB].resize(inst[blockType::CLB].size(), std::vector<int>(0));
    // RAM
    inst_to_net[blockType::RAM].resize(inst[blockType::RAM].size(), std::vector<int>(0));
    // DSP
    inst_to_net[blockType::DSP].resize(inst[blockType::DSP].size(), std::vector<int>(0));

    for (int net_id = 0; net_id < net_vec.size(); ++net_id) {
        for (auto &cell : net->inst_vec) {
            if (cell->type == blockType::IO) continue;
            else inst_to_net[cell->type][cell->id].emplace_back(net_id);
        }
    }
}

/*
void Solver::makeWindow() {
    // stable_sort(resource.begin(), resource.end(), [](const Block* lhs, const Block* rhs){
    //     return lhs->center_x < rhs->center_x;
    // });
    // stable_sort(resource.begin(), resource.end(), [](const Block* lhs, const Block* rhs){
    //     return lhs->center_y < rhs->center_y;
    // });
    std::vector<Block*> tmp = resource;
    stable_sort(tmp.begin(), tmp.end(), [](const Block* lhs, const Block* rhs){
        return lhs->center_y < rhs->center_y;
    });
    stable_sort(tmp.begin(), tmp.end(), [](const Block* lhs, const Block* rhs){
        return lhs->center_x < rhs->center_x;
    });

    int flag[4] = {2, 2, 2, 2};
    for (auto &it : tmp) {
        if (!flag[blockType::CLB] && !flag[blockType::RAM] && !flag[blockType::DSP]) break;
        if (!flag[it->type]) continue;
        if (it->type == blockType::CLB) {
            if (flag[blockType::CLB] == 2) {
                min_y[blockType::CLB] = it->center_y;
                height[blockType::CLB] -= it->center_y;
            } else if (flag[blockType::CLB] == 1) {
                height[blockType::CLB] += it->center_y;
            }
            flag[blockType::CLB]--;
        } else if (it->type == blockType::RAM) {
            if (flag[blockType::RAM] == 2) {
                min_y[blockType::RAM] = it->center_y;
                height[blockType::RAM] -= it->center_y;
            } else if (flag[blockType::RAM] == 1) {
                height[blockType::RAM] += it->center_y;
            }
            flag[blockType::RAM]--;
        } else if (it->type == blockType::DSP) {
            if (flag[blockType::DSP] == 2) {
                min_y[blockType::DSP] = it->center_y;
                height[blockType::DSP] -= it->center_y;
            } else if (flag[blockType::DSP] == 1) {
                height[blockType::DSP] += it->center_y;
            }
            flag[blockType::DSP]--;
        }
    }
    // std::cout << "CLB min y : " << min_y[blockType::CLB] << " RAM min y : " << min_y[blockType::RAM] << " DSP min y : " << min_y[blockType::DSP] << std::endl;
    // std::cout << "CLB height : " << height[blockType::CLB] << " RAM height : " << height[blockType::RAM] << " DSP height : " <<height[blockType::DSP] << std::endl;
    
    int max_h = std::max(height[blockType::CLB], std::max(height[blockType::RAM], height[blockType::DSP]));
    unit[blockType::CLB] = max_h / height[blockType::CLB];
    unit[blockType::RAM] = max_h / height[blockType::RAM];
    unit[blockType::DSP] = max_h / height[blockType::DSP];
    // std::cout << "CLB_unit: " << unit[blockType::CLB] << " RAM_unit: " << unit[blockType::RAM] << " DSP_unit: " << unit[blockType::DSP] << std::endl;
    
    // stable_sort(resource.begin(), resource.end(), [](const Block* lhs, const Block* rhs){
    //     int c1 = (lhs->center_y - min_y[lhs->type]) / height[lhs->type];
    //     int c2 = (rhs->center_y - min_y[rhs->type]) / height[rhs->type];
    //     if (c1 == c2) {
    //         return lhs->center_x < rhs->center_x;
    //     } else {
    //         return c1 < c2;
    //     }
    // });

    for (auto &it : resource) {
        it->layer = (it->center_y - min_y[it->type]) / height[it->type];
        it->unit = it->layer / unit[it->type];
    }

    stable_sort(resource.begin(), resource.end(), [](const Block* lhs, const Block* rhs){
        if (lhs->unit != rhs->unit) {
            return lhs->unit < rhs->unit;
        } else {
            if (lhs->layer != rhs->layer) return lhs->layer < rhs->layer;
            else return lhs->center_x < rhs->center_x;
        }
    });
    
    for (auto &it : resource) {
        if (it->type == blockType::CLB) std::cout << "[CLB]";
        else if (it->type == blockType::RAM) std::cout << "[RAM]";
        else if (it->type == blockType::DSP) std::cout << "[DSP]";
        else std::cout << "[ERROR] make window - 1\n";
        std::cout << " unit: " << it->unit << " layer: " << it->layer << " x: " << it->center_x << std::endl;

        // if (it->type == blockType::CLB) {
        //     int row_cnt = (it->center_y - min_y[blockType::CLB]) / height[blockType::CLB];
        //     if (row_cnt > CLB_row_cnt) CLB_row_cnt = row_cnt;
        // } else if (it->type == blockType::RAM) {
        //     int row_cnt = (it->center_y - min_y[blockType::RAM]) / height[blockType::RAM];
        //     if (row_cnt > RAM_row_cnt) RAM_row_cnt = row_cnt;
        // } else if (it->type == blockType::DSP) {
        //     int row_cnt = (it->center_y - min_y[blockType::DSP]) / height[blockType::DSP];
        //     if (row_cnt > DSP_row_cnt) DSP_row_cnt = row_cnt;
        // } else {
        //     std::cout << "[ERROR] make window - 1\n";
        // }
    }   

    std::vector<Block*> single_window[4];
    int u = resource[0]->unit;
    for (auto &it : resource) {
        if (it->unit != u) {
            window[blockType::CLB].emplace_back(single_window[blockType::CLB]);
            window[blockType::RAM].emplace_back(single_window[blockType::RAM]);
            window[blockType::DSP].emplace_back(single_window[blockType::DSP]);
            single_window[blockType::CLB].clear();
            single_window[blockType::RAM].clear();
            single_window[blockType::DSP].clear();
            u = it->unit;
            single_window[it->type].emplace_back(it);
        } else {
            single_window[it->type].emplace_back(it);
        }
    }
    // std::cout << "CLB: \n";
    // std::cout << window[blockType::CLB].size() << std::endl;
    // for (auto &it : window[blockType::CLB]) std::cout << it.size() << ' ';
    // std::cout << std::endl;
    // std::cout << "RAM: \n";
    // std::cout << window[blockType::RAM].size() << std::endl;
    // for (auto &it : window[blockType::RAM]) std::cout << it.size() << ' ';
    // std::cout << std::endl;
    // std::cout << "DSP: \n";
    // std::cout << window[blockType::DSP].size() << std::endl;
    // for (auto &it : window[blockType::DSP]) std::cout << it.size() << ' ';
    // std::cout << std::endl;
}
*/

void Solver::fitness(gene& g) {
    double min_x, min_y, max_x, max_y;
    double fitness = 0;
    for (auto &net : net_vec) {
        min_x = INT_MAX;
        min_y = INT_MAX;
        max_x = INT_MIN;
        max_y = INT_MIN;
        for (auto &inst : net->inst_vec) {
            if (inst->type == blockType::IO) {
                if (inst->center_x < min_x) min_x = inst->center_x;
                if (inst->center_x > max_x) max_x = inst->center_x;
                if (inst->center_y < min_y) min_y = inst->center_y;
                if (inst->center_y > max_y) max_y = inst->center_y;
            } else {
                if (resource[inst->type][g.resource_permu[inst->type][inst->id]]->center_x < min_x) min_x = resource[inst->type][g.resource_permu[inst->type][inst->id]]->center_x;
                if (resource[inst->type][g.resource_permu[inst->type][inst->id]]->center_x > max_x) max_x = resource[inst->type][g.resource_permu[inst->type][inst->id]]->center_x;
                if (resource[inst->type][g.resource_permu[inst->type][inst->id]]->center_y < min_y) min_y = resource[inst->type][g.resource_permu[inst->type][inst->id]]->center_y;
                if (resource[inst->type][g.resource_permu[inst->type][inst->id]]->center_y > max_y) max_y = resource[inst->type][g.resource_permu[inst->type][inst->id]]->center_y;
            }
        }
        fitness += (max_x - min_x) + (max_y - min_y);
    }
    g.fitness = fitness;
}

void Solver::init_pop() {
    for (size_t i = 0; i < POP_SIZE; ++i) {
        gene p;
        // CLB
        for (size_t j = 0; j < resource[blockType::CLB].size(); ++j) {
            p.resource_permu[blockType::CLB].emplace_back(j);
        }
        for(size_t j = resource[blockType::CLB].size() - 1; j >= 1; --j) {
            std::swap(p.resource_permu[blockType::CLB][j], p.resource_permu[blockType::CLB][rand() % j]);
        }
        // RAM
        for (size_t j = 0; j < resource[blockType::RAM].size(); ++j) {
            p.resource_permu[blockType::RAM].emplace_back(j);
        }
        for(size_t j = resource[blockType::RAM].size() - 1; j >= 1; --j) {
            std::swap(p.resource_permu[blockType::RAM][j], p.resource_permu[blockType::RAM][rand() % j]);
        }
        // DSP
        for (size_t j = 0; j < resource[blockType::DSP].size(); ++j) {
            p.resource_permu[blockType::DSP].emplace_back(j);
        }
        for(size_t j = resource[blockType::DSP].size() - 1; j >= 1; --j) {
            std::swap(p.resource_permu[blockType::DSP][j], p.resource_permu[blockType::DSP][rand() % j]);
        }

        pool.emplace_back(p);
    }

    // for (size_t i = 0; i < pool.size(); ++i) {
    //     std::cout << "CLB: " << std::endl;
    //     for (size_t j = 0; j < pool[i].resource_permu[blockType::CLB].size(); ++j) {
    //         std::cout << pool[i].resource_permu[blockType::CLB][j] << ' ';
    //     }
    //     std::cout << std::endl;
    //     std::cout << "RAM: " << std::endl;
    //     for (size_t j = 0; j < pool[i].resource_permu[blockType::RAM].size(); ++j) {
    //         std::cout << pool[i].resource_permu[blockType::RAM][j] << ' ';
    //     }
    //     std::cout << std::endl;
    //     std::cout << "DSP: " << std::endl;
    //     for (size_t j = 0; j < pool[i].resource_permu[blockType::DSP].size(); ++j) {
    //         std::cout << pool[i].resource_permu[blockType::DSP][j] << ' ';
    //     }
    //     std::cout << std::endl;
    // }

    for (auto &g : pool) fitness(g);
    std::sort(pool.begin(), pool.end(), [](const gene& lhs, const gene& rhs){
        return lhs.fitness < rhs.fitness;
    });

    std::cout << pool[0].fitness << std:: endl;
}

void Solver::parent_selection(parents& parent) {
    std::sort(pool.begin(), pool.end(), [](const gene& lhs, const gene& rhs){
        return lhs.fitness < rhs.fitness;
    });
    int idx1, idx2;
    idx1 = idx2 = INT_MAX;
    for (int i = 0, tmp1, tmp2; i < K; ++i) {
        tmp1 = rand() % pool.size();
        tmp2 = rand() % pool.size();
        idx1 = std::min(idx1, tmp1);
        idx2 = std::min(idx2, tmp2);
    }
    parent.first = pool[idx1];
    parent.second = pool[(idx1 != idx2) ? idx2 : idx2 + 1];
}

void Solver::crossover(parents& parent, std::vector<gene>& offspring, int type) {
    gene child1, child2;
    child1.resource_permu[blockType::CLB].resize(resource[blockType::CLB].size(), 0);
    child1.resource_permu[blockType::RAM].resize(resource[blockType::RAM].size(), 0);
    child1.resource_permu[blockType::DSP].resize(resource[blockType::DSP].size(), 0);
    child2.resource_permu[blockType::CLB].resize(resource[blockType::CLB].size(), 0);
    child2.resource_permu[blockType::RAM].resize(resource[blockType::RAM].size(), 0);
    child2.resource_permu[blockType::DSP].resize(resource[blockType::DSP].size(), 0);
    
    if (type == 0) { // order crossover
        int tmp1, tmp2;
        // CLB
        child1.resource_permu[blockType::CLB] = parent.first.resource_permu[blockType::CLB];
        child2.resource_permu[blockType::CLB] = parent.second.resource_permu[blockType::CLB];
        do {
            tmp1 = rand() % resource[blockType::CLB].size();
            tmp2 = rand() % resource[blockType::CLB].size();
        } while (tmp1 == tmp2);
        if (tmp1 > tmp2) std::swap(tmp1, tmp2);
        st1.clear();
        st2.clear();
        for (size_t i = tmp1; i <= tmp2; ++i) {
            st1.insert(parent.first.resource_permu[blockType::CLB][i]);
            st2.insert(parent.second.resource_permu[blockType::CLB][i]);
        }
        int idx2 = tmp2, idx1 = tmp2;
        for (int i = tmp2 + 1; i < resource[blockType::CLB].size(); ++i) {
            idx2++;
            while (st1.find(parent.second.resource_permu[blockType::CLB][mod_table[blockType::CLB][idx2]]) != st1.end()) {
                idx2++;
            }
            child1.resource_permu[blockType::CLB][i] = parent.second.resource_permu[blockType::CLB][mod_table[blockType::CLB][idx2]];

            idx1++;
            while (st2.find(parent.first.resource_permu[blockType::CLB][mod_table[blockType::CLB][idx1]]) != st2.end()) {
                idx1++;
            }
            child2.resource_permu[blockType::CLB][i] = parent.first.resource_permu[blockType::CLB][mod_table[blockType::CLB][idx1]];
        }
        for (int i = 0; i < tmp1; ++i) {
            idx2++;
            while (st1.find(parent.second.resource_permu[blockType::CLB][mod_table[blockType::CLB][idx2]]) != st1.end()) {
                idx2++;
            }
            child1.resource_permu[blockType::CLB][i] = parent.second.resource_permu[blockType::CLB][mod_table[blockType::CLB][idx2]];

            idx1++;
            while (st2.find(parent.first.resource_permu[blockType::CLB][mod_table[blockType::CLB][idx1]]) != st2.end()) {
                idx1++;
            }
            child2.resource_permu[blockType::CLB][i] = parent.first.resource_permu[blockType::CLB][mod_table[blockType::CLB][idx1]];
        }

        // RAM
        child1.resource_permu[blockType::RAM] = parent.first.resource_permu[blockType::RAM];
        child2.resource_permu[blockType::RAM] = parent.second.resource_permu[blockType::RAM];
        do {
            tmp1 = rand() % resource[blockType::RAM].size();
            tmp2 = rand() % resource[blockType::RAM].size();
        } while (tmp1 == tmp2);
        if (tmp1 > tmp2) std::swap(tmp1, tmp2);
        st1.clear();
        st2.clear();
        for (size_t i = tmp1; i <= tmp2; ++i) {
            st1.insert(parent.first.resource_permu[blockType::RAM][i]);
            st2.insert(parent.second.resource_permu[blockType::RAM][i]);
        }
        idx2 = tmp2, idx1 = tmp2;
        for (int i = tmp2 + 1; i < resource[blockType::RAM].size(); ++i) {
            idx2++;
            while (st1.find(parent.second.resource_permu[blockType::RAM][mod_table[blockType::RAM][idx2]]) != st1.end()) {
                idx2++;
            }
            child1.resource_permu[blockType::RAM][i] = parent.second.resource_permu[blockType::RAM][mod_table[blockType::RAM][idx2]];

            idx1++;
            while (st2.find(parent.first.resource_permu[blockType::RAM][mod_table[blockType::RAM][idx1]]) != st2.end()) {
                idx1++;
            }
            child2.resource_permu[blockType::RAM][i] = parent.first.resource_permu[blockType::RAM][mod_table[blockType::RAM][idx1]];
        }
        for (int i = 0; i < tmp1; ++i) {
            idx2++;
            while (st1.find(parent.second.resource_permu[blockType::RAM][mod_table[blockType::RAM][idx2]]) != st1.end()) {
                idx2++;
            }
            child1.resource_permu[blockType::RAM][i] = parent.second.resource_permu[blockType::RAM][mod_table[blockType::RAM][idx2]];

            idx1++;
            while (st2.find(parent.first.resource_permu[blockType::RAM][mod_table[blockType::RAM][idx1]]) != st2.end()) {
                idx1++;
            }
            child2.resource_permu[blockType::RAM][i] = parent.first.resource_permu[blockType::RAM][mod_table[blockType::RAM][idx1]];
        }

        // DSP
        child1.resource_permu[blockType::DSP] = parent.first.resource_permu[blockType::DSP];
        child2.resource_permu[blockType::DSP] = parent.second.resource_permu[blockType::DSP];
        do {
            tmp1 = rand() % resource[blockType::DSP].size();
            tmp2 = rand() % resource[blockType::DSP].size();
        } while (tmp1 == tmp2);
        if (tmp1 > tmp2) std::swap(tmp1, tmp2);
        st1.clear();
        st2.clear();
        for (size_t i = tmp1; i <= tmp2; ++i) {
            st1.insert(parent.first.resource_permu[blockType::DSP][i]);
            st2.insert(parent.second.resource_permu[blockType::DSP][i]);
        }
        idx2 = tmp2, idx1 = tmp2;
        for (int i = tmp2 + 1; i < resource[blockType::DSP].size(); ++i) {
            idx2++;
            while (st1.find(parent.second.resource_permu[blockType::DSP][mod_table[blockType::DSP][idx2]]) != st1.end()) {
                idx2++;
            }
            child1.resource_permu[blockType::DSP][i] = parent.second.resource_permu[blockType::DSP][mod_table[blockType::DSP][idx2]];

            idx1++;
            while (st2.find(parent.first.resource_permu[blockType::DSP][mod_table[blockType::DSP][idx1]]) != st2.end()) {
                idx1++;
            }
            child2.resource_permu[blockType::DSP][i] = parent.first.resource_permu[blockType::DSP][mod_table[blockType::DSP][idx1]];
        }
        for (int i = 0; i < tmp1; ++i) {
            idx2++;
            while (st1.find(parent.second.resource_permu[blockType::DSP][mod_table[blockType::DSP][idx2]]) != st1.end()) {
                idx2++;
            }
            child1.resource_permu[blockType::DSP][i] = parent.second.resource_permu[blockType::DSP][mod_table[blockType::DSP][idx2]];

            idx1++;
            while (st2.find(parent.first.resource_permu[blockType::DSP][mod_table[blockType::DSP][idx1]]) != st2.end()) {
                idx1++;
            }
            child2.resource_permu[blockType::DSP][i] = parent.first.resource_permu[blockType::DSP][mod_table[blockType::DSP][idx1]];
        }
    } else if (type == 1) { // cycle crossover
        // CLB
        child1.resource_permu[blockType::CLB] = parent.first.resource_permu[blockType::CLB];
        child2.resource_permu[blockType::CLB] = parent.second.resource_permu[blockType::CLB];
        int turn = 0, val1;
        m.clear();
        for (int i = 0; i < resource[blockType::CLB].size(); ++i) {
            if (!m.count(child1.resource_permu[blockType::CLB][i])) m[child1.resource_permu[blockType::CLB][i]] = i;
        }
        while (!m.empty()) {
            val1 = m.begin()->first; // child1's value = child2's index
            while (1) {
                if (!m.count(val1)) break; // a cycle has been formed
                if (turn & 1) std::swap(child1.resource_permu[blockType::CLB][m[val1]], child2.resource_permu[blockType::CLB][val1]);
                val1 = child1.resource_permu[blockType::CLB][m[val1]];
                m.erase(val1);
            }
            turn++;
        }
        
        // RAM
        child1.resource_permu[blockType::RAM] = parent.first.resource_permu[blockType::RAM];
        child2.resource_permu[blockType::RAM] = parent.second.resource_permu[blockType::RAM];
        turn = 0;
        m.clear();
        for (int i = 0; i < resource[blockType::RAM].size(); ++i) {
            if (!m.count(child1.resource_permu[blockType::RAM][i])) m[child1.resource_permu[blockType::RAM][i]] = i;
        }
        while (!m.empty()) {
            val1 = m.begin()->first; // child1's value = child2's index
            while (1) {
                if (!m.count(val1)) break; // a cycle has been formed
                if (turn & 1) std::swap(child1.resource_permu[blockType::RAM][m[val1]], child2.resource_permu[blockType::RAM][val1]);
                val1 = child1.resource_permu[blockType::RAM][m[val1]];
                m.erase(val1);
            }
            turn++;
        }

        // DSP
        child1.resource_permu[blockType::DSP] = parent.first.resource_permu[blockType::DSP];
        child2.resource_permu[blockType::DSP] = parent.second.resource_permu[blockType::DSP];
        turn = 0;
        m.clear();
        for (int i = 0; i < resource[blockType::DSP].size(); ++i) {
            if (!m.count(child1.resource_permu[blockType::DSP][i])) m[child1.resource_permu[blockType::DSP][i]] = i;
        }
        while (!m.empty()) {
            val1 = m.begin()->first; // child1's value = child2's index
            while (1) {
                if (!m.count(val1)) break; // a cycle has been formed
                if (turn & 1) std::swap(child1.resource_permu[blockType::DSP][m[val1]], child2.resource_permu[blockType::DSP][val1]);
                val1 = child1.resource_permu[blockType::DSP][m[val1]];
                m.erase(val1);
            }
            turn++;
        }
    }
    offspring.emplace_back(child1);
    offspring.emplace_back(child2);
}

void Solver::mutation(std::vector<gene>& offspring, int type) {
    if (type == 0) { // swap mutation
        int tmp1, tmp2;
        // CLB
        tmp1 = rand() % inst[blockType::CLB].size();
        tmp2 = rand() % inst[blockType::CLB].size();
        std::swap(offspring[offspring.size() - 1].resource_permu[blockType::CLB][tmp1], offspring[offspring.size() - 1].resource_permu[blockType::CLB][tmp2]);
        tmp1 = rand() % inst[blockType::CLB].size();
        tmp2 = rand() % inst[blockType::CLB].size();
        std::swap(offspring[offspring.size() - 2].resource_permu[blockType::CLB][tmp1], offspring[offspring.size() - 2].resource_permu[blockType::CLB][tmp2]);

        // RAM
        tmp1 = rand() % inst[blockType::RAM].size();
        tmp2 = rand() % inst[blockType::RAM].size();
        std::swap(offspring[offspring.size() - 1].resource_permu[blockType::RAM][tmp1], offspring[offspring.size() - 1].resource_permu[blockType::RAM][tmp2]);
        tmp1 = rand() % inst[blockType::RAM].size();
        tmp2 = rand() % inst[blockType::RAM].size();
        std::swap(offspring[offspring.size() - 2].resource_permu[blockType::RAM][tmp1], offspring[offspring.size() - 2].resource_permu[blockType::RAM][tmp2]);

        // DSP
        tmp1 = rand() % inst[blockType::DSP].size();
        tmp2 = rand() % inst[blockType::DSP].size();
        std::swap(offspring[offspring.size() - 1].resource_permu[blockType::DSP][tmp1], offspring[offspring.size() - 1].resource_permu[blockType::DSP][tmp2]);
        tmp1 = rand() % inst[blockType::DSP].size();
        tmp2 = rand() % inst[blockType::DSP].size();
        std::swap(offspring[offspring.size() - 2].resource_permu[blockType::DSP][tmp1], offspring[offspring.size() - 2].resource_permu[blockType::DSP][tmp2]);

        fitness(offspring[offspring.size() - 1]);
        fitness(offspring[offspring.size() - 2]);
    }
}

void Solver::survivor_selection(std::vector<gene>& new_genes) {
    pool.insert(pool.end(), new_genes.begin(), new_genes.end());
    std::sort(pool.begin(), pool.end(), [](const gene& lhs, const gene& rhs){
        return lhs.fitness < rhs.fitness;
    });
    pool = std::vector<gene>(pool.begin(), pool.begin() + POP_SIZE);
}

void Solver::genetic_algorithm() {
    parents parent;
    std::vector<gene> offspring;
    for (int i = 0; i < TERMINATION; ++i) {
        for (int j = 0; j < POP_SIZE / 2; ++j) {
            parent_selection(parent);
            // std::cout << "finish parent selection\n";
            crossover(parent, offspring, 1);
            // std::cout << "finish crossover\n";
            mutation(offspring, 0);
            // std::cout << "finish mutation\n";
        }
        survivor_selection(offspring);
        offspring.clear();
        std::cout << "[round " << i+1 << "]: " << pool[0].fitness << std::endl;
    }
}

void Solver::output_file(char* output_file) {
    std::ofstream fout;
    fout.open(output_file);
    for (size_t i = 0; i < inst[blockType::CLB].size(); ++i) {
        fout << inst[blockType::CLB][i]->name << " " << resource[blockType::CLB][pool[0].resource_permu[blockType::CLB][i]]->name << std::endl;
    }
    for (size_t i = 0; i < inst[blockType::RAM].size(); ++i) {
        fout << inst[blockType::RAM][i]->name << " " << resource[blockType::RAM][pool[0].resource_permu[blockType::RAM][i]]->name << std::endl;
    }
    for (size_t i = 0; i < inst[blockType::DSP].size(); ++i) {
        fout << inst[blockType::DSP][i]->name << " " << resource[blockType::DSP][pool[0].resource_permu[blockType::DSP][i]]->name << std::endl;
    }
    fout.close();
}


