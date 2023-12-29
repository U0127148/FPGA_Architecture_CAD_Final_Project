#include "solver.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <random>
#include <climits>
#include <unordered_set>
#include <unordered_map>
#include <cmath>
const int POP_SIZE = 100;
const int K = 2;
double min_y[4];
double height[4];
std::unordered_map<int, int> m; // (child1 value, child1 index)
std::vector<std::vector<Block*>> resource_grid[4];

void Solver::read(char *argv[]) {
    std::string line;

    int idx[4] = {0, 0, 0, 0};
    // first input file (resource)
    std::ifstream fin;
    fin.open(argv[1]);
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
                    std::cout << "[ERROR] read FIRST input file\n"; 
                }
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
                resource_block->id = id;
            }
            n++;
        }
    }
    fin.close();

    // for (auto& cell : resource[blockType::CLB]) std::cout << cell->id << ' ';
    // std::cout << std::endl;
    // for (auto& cell : resource[blockType::RAM]) std::cout << cell->id << ' ';
    // std::cout << std::endl;
    // for (auto& cell : resource[blockType::DSP]) std::cout << cell->id << ' ';
    // std::cout << std::endl;

    idx[blockType::CLB] = 0;
    idx[blockType::RAM] = 0;
    idx[blockType::DSP] = 0;
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
}

void Solver::create_grid() {
    std::vector<Block*> tmp[4];
    tmp[blockType::CLB] = resource[blockType::CLB];
    tmp[blockType::RAM] = resource[blockType::RAM];
    tmp[blockType::DSP] = resource[blockType::DSP];
    // sort the resources of CLB, RAM, DSP by x-coordinate then by y-coordinate
    stable_sort(tmp[blockType::CLB].begin(), tmp[blockType::CLB].end(), [](const Block* lhs, const Block* rhs){
        return lhs->center_y < rhs->center_y;
    });
    stable_sort(tmp[blockType::RAM].begin(), tmp[blockType::RAM].end(), [](const Block* lhs, const Block* rhs){
        return lhs->center_y < rhs->center_y;
    });
    stable_sort(tmp[blockType::DSP].begin(), tmp[blockType::DSP].end(), [](const Block* lhs, const Block* rhs){
        return lhs->center_y < rhs->center_y;
    });
    // for (auto &it : tmp[blockType::DSP]) std::cout << it->center_x << ' ' << it->center_y << std::endl;

    int flag = 2;
    // get the min y-coordinate and the height of CLB, RAM, DSP
    for (auto &it : tmp[blockType::CLB]) {
        if (flag == 2) {
            min_y[blockType::CLB] = it->center_y;
            flag--;
        } else if (flag == 1) {
            if (it->center_y != min_y[blockType::CLB]) {
                height[blockType::CLB] = it->center_y - min_y[blockType::CLB];
                flag--;
            }
        } else {
            break;
        }
    }
    flag = 2;
    for (auto &it : tmp[blockType::RAM]) {
        if (flag == 2) {
            min_y[blockType::RAM] = it->center_y;
            flag--;
        } else if (flag == 1) {
            if (it->center_y != min_y[blockType::RAM]) {
                height[blockType::RAM] = it->center_y - min_y[blockType::RAM];
                flag--;
            }
        } else {
            break;
        }
    }
    flag = 2;
    for (auto &it : tmp[blockType::DSP]) {
        if (flag == 2) {
            min_y[blockType::DSP] = it->center_y;
            flag--;
        } else if (flag == 1) {
            if (it->center_y != min_y[blockType::DSP]) {
                height[blockType::DSP] = it->center_y - min_y[blockType::DSP];
                flag--;
            }
        } else {
            break;
        }
    }
    // std::cout << "CLB min y : " << min_y[blockType::CLB] << " RAM min y : " << min_y[blockType::RAM] << " DSP min y : " << min_y[blockType::DSP] << std::endl;
    // std::cout << "CLB height : " << height[blockType::CLB] << " RAM height : " << height[blockType::RAM] << " DSP height : " <<height[blockType::DSP] << std::endl;

    // put the block in the same height for every tmp type
    std::vector<Block*> v;
    double h = min_y[blockType::CLB];
    for (auto &it : tmp[blockType::CLB]) {
        if (it->center_y != h) {
            // std::cout << " new height " << it->center_x << ' ' << it->center_y << std::endl;
            stable_sort(v.begin(), v.end(), [](const Block* lhs, const Block* rhs){
                return lhs->center_x < rhs->center_x;
            });
            resource_grid[blockType::CLB].emplace_back(v);
            v.clear();
            v.emplace_back(it);
            h = it->center_y;
        } else {
            v.emplace_back(it);
        }
    }
    stable_sort(v.begin(), v.end(), [](const Block* lhs, const Block* rhs){
        return lhs->center_x < rhs->center_x;
    });
    resource_grid[blockType::CLB].emplace_back(v);

    v.clear();
    h = tmp[blockType::RAM][0]->center_y;
    for (auto &it : tmp[blockType::RAM]) {
        if (it->center_y != h) {
            stable_sort(v.begin(), v.end(), [](const Block* lhs, const Block* rhs){
                return lhs->center_x < rhs->center_x;
            });
            resource_grid[blockType::RAM].emplace_back(v);
            v.clear();
            v.emplace_back(it);
            h = it->center_y;
        } else {
            v.emplace_back(it);
        }
    }
    stable_sort(v.begin(), v.end(), [](const Block* lhs, const Block* rhs){
        return lhs->center_x < rhs->center_x;
    });
    resource_grid[blockType::RAM].emplace_back(v);

    v.clear();
    h = tmp[blockType::DSP][0]->center_y;
    for (auto &it : tmp[blockType::DSP]) {
        if (it->center_y != h) {
            stable_sort(v.begin(), v.end(), [](const Block* lhs, const Block* rhs){
                return lhs->center_x < rhs->center_x;
            });
            resource_grid[blockType::DSP].emplace_back(v);
            v.clear();
            v.emplace_back(it);
            h = it->center_y;
        } else {
            v.emplace_back(it);
        }
    }
    stable_sort(v.begin(), v.end(), [](const Block* lhs, const Block* rhs){
        return lhs->center_x < rhs->center_x;
    });
    resource_grid[blockType::DSP].emplace_back(v);
}

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

int Solver::find_available_block(int type, int r, int c, int p) {
    int dist = 1;
    int dx[8] = {0, 1, 0, -1, 1, 1, -1, -1};
    int dy[8] = {1, 0, -1, 0, 1, -1, -1, 1};
    while (1) {
        for (size_t i = 0; i < 8; ++i) {
            if ((r + dist * dy[i] >= resource_grid[type].size()) || (r + dist * dy[i] < 0) || (c + dist * dx[i] >= resource_grid[type][r].size()) || (c + dist * dx[i] < 0)) continue;
            if (resource_grid[type][r + dist * dy[i]][c + dist * dx[i]]->sel == p) {
                resource_grid[type][r + dist * dy[i]][c + dist * dx[i]]->sel++;
                return resource_grid[type][r + dist * dy[i]][c + dist * dx[i]]->id;
            }
        }
        dist++;
    }
    return 0;
}

int Solver::find_block(int type, int inst_idx, int p) {
    if (type == blockType::CLB) {
        int r = (inst[blockType::CLB][inst_idx]->center_y < min_y[blockType::CLB]) ? 0 : (inst[blockType::CLB][inst_idx]->center_y - min_y[blockType::CLB]) / height[blockType::CLB];
        if (r >= resource_grid[blockType::CLB].size() - 1) r = resource_grid[blockType::CLB].size() - 1;
        else r = (abs(resource_grid[blockType::CLB][r][0]->center_y - inst[blockType::CLB][inst_idx]->center_y) < abs(resource_grid[blockType::CLB][r + 1][0]->center_y - inst[blockType::CLB][inst_idx]->center_y)) ? r : r + 1;
        int c;
        if (inst[blockType::CLB][inst_idx]->center_x < resource_grid[blockType::CLB][r][0]->center_x) c = 0;
        else {
            for (size_t i = 0; i < resource_grid[blockType::CLB][r].size(); ++i) {
                if (resource_grid[blockType::CLB][r][i]->center_x > inst[blockType::CLB][inst_idx]->center_x) {
                    c = (abs(resource_grid[blockType::CLB][r][i]->center_x - inst[blockType::CLB][inst_idx]->center_x) < abs(resource_grid[blockType::CLB][r][i - 1]->center_x - inst[blockType::CLB][inst_idx]->center_x)) ? i : i - 1;
                    break;
                }
            }
        }
        if (resource_grid[blockType::CLB][r][c]->sel == p) {
            resource_grid[blockType::CLB][r][c]->sel++;
            return resource_grid[blockType::CLB][r][c]->id;
        } else {
            return find_available_block(type, r, c, p);
        }
    } else if (type == blockType::RAM) {
        int r = (inst[blockType::RAM][inst_idx]->center_y < min_y[blockType::RAM]) ? 0 : (inst[blockType::RAM][inst_idx]->center_y - min_y[blockType::RAM]) / height[blockType::RAM];
        if (r >= resource_grid[blockType::RAM].size() - 1) r = resource_grid[blockType::RAM].size() - 1;
        else r = (abs(resource_grid[blockType::RAM][r][0]->center_y - inst[blockType::RAM][inst_idx]->center_y) < abs(resource_grid[blockType::RAM][r + 1][0]->center_y - inst[blockType::RAM][inst_idx]->center_y)) ? r : r + 1;
        int c;
        if (inst[blockType::RAM][inst_idx]->center_x < resource_grid[blockType::RAM][r][0]->center_x) c = 0;
        else {
            for (size_t i = 0; i < resource_grid[blockType::RAM][r].size(); ++i) {
                if (resource_grid[blockType::RAM][r][i]->center_x > inst[blockType::RAM][inst_idx]->center_x) {
                    c = (abs(resource_grid[blockType::RAM][r][i]->center_x - inst[blockType::RAM][inst_idx]->center_x) < abs(resource_grid[blockType::RAM][r][i - 1]->center_x - inst[blockType::RAM][inst_idx]->center_x)) ? i : i - 1;
                    break;
                }
            }
        }
        if (resource_grid[blockType::RAM][r][c]->sel == p) {
            resource_grid[blockType::RAM][r][c]->sel++;
            return resource_grid[blockType::RAM][r][c]->id;
        } else {
            return find_available_block(type, r, c, p);
        }
    } else {
        int r = (inst[blockType::DSP][inst_idx]->center_y < min_y[blockType::DSP]) ? 0 : (inst[blockType::DSP][inst_idx]->center_y - min_y[blockType::DSP]) / height[blockType::DSP];
        if (r >= resource_grid[blockType::DSP].size() - 1) r = resource_grid[blockType::DSP].size() - 1;
        else r = (abs(resource_grid[blockType::DSP][r][0]->center_y - inst[blockType::DSP][inst_idx]->center_y) < abs(resource_grid[blockType::DSP][r + 1][0]->center_y - inst[blockType::DSP][inst_idx]->center_y)) ? r : r + 1;
        int c;
        if (inst[blockType::DSP][inst_idx]->center_x < resource_grid[blockType::DSP][r][0]->center_x) c = 0;
        else {
            for (size_t i = 0; i < resource_grid[blockType::DSP][r].size(); ++i) {
                if (resource_grid[blockType::DSP][r][i]->center_x > inst[blockType::DSP][inst_idx]->center_x) {
                    c = (abs(resource_grid[blockType::DSP][r][i]->center_x - inst[blockType::DSP][inst_idx]->center_x) < abs(resource_grid[blockType::DSP][r][i - 1]->center_x - inst[blockType::DSP][inst_idx]->center_x)) ? i : i - 1;
                    break;
                }
            }
        }
        if (resource_grid[blockType::DSP][r][c]->sel == p) {
            resource_grid[blockType::DSP][r][c]->sel++;
            return resource_grid[blockType::DSP][r][c]->id;
        } else {
            return find_available_block(type, r, c, p);
        }
    }
    return 0;
}

void Solver::init_pop_from_GP() {
    std::unordered_set<int> tmp;
    std::vector<int> v1, v2, v3;
    for (int i = 0; i < inst[blockType::CLB].size(); ++i) v1.emplace_back(i);
    for (int i = 0; i < inst[blockType::RAM].size(); ++i) v2.emplace_back(i);
    for (int i = 0; i < inst[blockType::DSP].size(); ++i) v3.emplace_back(i);

    int remain_idx;
    std::vector<gene> tmp_pool;
    for (size_t i = 0; i < POP_SIZE + 50; ++i) {
        gene p;
        // CLB
        p.resource_permu[blockType::CLB].resize(resource[blockType::CLB].size());
        for (size_t j = inst[blockType::CLB].size() - 1; j >= 1; --j) {
            std::swap(v1[j], v1[rand() % j]);
        }
        for (int j = 0, idx; j < inst[blockType::CLB].size(); ++j) {
            idx = find_block(blockType::CLB, v1[j], i);
            p.resource_permu[blockType::CLB][v1[j]] = idx;
            tmp.insert(idx);
        }
        remain_idx = 0;
        for (size_t j = inst[blockType::CLB].size(); j < resource[blockType::CLB].size(); ++j) {
            while (tmp.find(remain_idx) != tmp.end()) remain_idx++;
            resource[blockType::CLB][remain_idx]->sel++;
            p.resource_permu[blockType::CLB][j] = remain_idx++;
        }
        for (size_t j = resource[blockType::CLB].size() - 1; j >= inst[blockType::CLB].size(); --j) {
            int random_tmp = rand() % (j - inst[blockType::CLB].size() + 1) + inst[blockType::CLB].size();
            std::swap(p.resource_permu[blockType::CLB][j], p.resource_permu[blockType::CLB][random_tmp]);
        }
        tmp.clear();

        // RAM
        p.resource_permu[blockType::RAM].resize(resource[blockType::RAM].size());
        for (size_t j = inst[blockType::RAM].size() - 1; j >= 1; --j) {
            std::swap(v2[j], v2[rand() % j]);
        }
        for (int j = 0, idx; j < inst[blockType::RAM].size(); ++j) {
            idx = find_block(blockType::RAM, v2[j], i);
            p.resource_permu[blockType::RAM][v2[j]] = idx;
            tmp.insert(idx);
        }
        remain_idx = 0;
        for (size_t j = inst[blockType::RAM].size(); j < resource[blockType::RAM].size(); ++j) {
            while (tmp.find(remain_idx) != tmp.end()) remain_idx++;
            resource[blockType::RAM][remain_idx]->sel++;
            p.resource_permu[blockType::RAM][j] = remain_idx++;
        }
        for (size_t j = resource[blockType::RAM].size() - 1; j >= inst[blockType::RAM].size(); --j) {
            int random_tmp = rand() % (j - inst[blockType::RAM].size() + 1) + inst[blockType::RAM].size();
            std::swap(p.resource_permu[blockType::RAM][j], p.resource_permu[blockType::RAM][random_tmp]);
        }
        tmp.clear();
        
        // DSP
        p.resource_permu[blockType::DSP].resize(resource[blockType::DSP].size());
        for (size_t j = inst[blockType::DSP].size() - 1; j >= 1; --j) {
            std::swap(v3[j], v3[rand() % j]);
        }
        for (int j = 0, idx; j < inst[blockType::DSP].size(); ++j) {
            idx = find_block(blockType::DSP, v3[j], i);
            p.resource_permu[blockType::DSP][v3[j]] = idx;
            tmp.insert(idx);
        }
        remain_idx = 0;
        for (size_t j = inst[blockType::DSP].size(); j < resource[blockType::DSP].size(); ++j) {
            while (tmp.find(remain_idx) != tmp.end()) remain_idx++;
            resource[blockType::DSP][remain_idx]->sel++;
            p.resource_permu[blockType::DSP][j] = remain_idx++;
        }
        for (size_t j = resource[blockType::DSP].size() - 1; j >= inst[blockType::DSP].size(); --j) {
            int random_tmp = rand() % (j - inst[blockType::DSP].size() + 1) + inst[blockType::DSP].size();
            std::swap(p.resource_permu[blockType::DSP][j], p.resource_permu[blockType::DSP][random_tmp]);
        }
        tmp.clear();

        tmp_pool.emplace_back(p);
    }

    for (auto &g : tmp_pool) fitness(g);
    std::sort(tmp_pool.begin(), tmp_pool.end(), [](const gene& lhs, const gene& rhs){
        return lhs.fitness < rhs.fitness;
    });
    pool.insert(pool.end(), tmp_pool.begin(), tmp_pool.begin() + POP_SIZE);

    std::cout << "Initial HPWL: " << pool[0].fitness << std:: endl;
}

void Solver::parent_selection(parents& parent) {
    std::sort(pool.begin(), pool.end(), [](const gene& lhs, const gene& rhs){
        return lhs.fitness < rhs.fitness;
    });
    int idx1, idx2;
    idx1 = idx2 = INT_MAX;
    for (int i = 0, tmp1, tmp2; i < K; ++i) {
        tmp1 = rand() % (POP_SIZE - 1);
        tmp2 = rand() % (POP_SIZE - 1);
        idx1 = std::min(idx1, tmp1);
        idx2 = std::min(idx2, tmp2);
    }
    parent.first = pool[idx1];
    parent.second = pool[(idx1 != idx2) ? idx2 : idx2 + 1];
}

void Solver::crossover(parents& parent, std::vector<gene>& offspring) {
    gene child1, child2;
    child1.resource_permu[blockType::CLB].resize(resource[blockType::CLB].size(), 0);
    child1.resource_permu[blockType::RAM].resize(resource[blockType::RAM].size(), 0);
    child1.resource_permu[blockType::DSP].resize(resource[blockType::DSP].size(), 0);
    child2.resource_permu[blockType::CLB].resize(resource[blockType::CLB].size(), 0);
    child2.resource_permu[blockType::RAM].resize(resource[blockType::RAM].size(), 0);
    child2.resource_permu[blockType::DSP].resize(resource[blockType::DSP].size(), 0);
    
    // cycle crossover
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
    offspring.emplace_back(child1);
    offspring.emplace_back(child2);
}

void Solver::mutation(std::vector<gene>& offspring) {
    // swap mutation
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

void Solver::survivor_selection(std::vector<gene>& new_genes) {
    pool.insert(pool.end(), new_genes.begin(), new_genes.end());
    std::sort(pool.begin(), pool.end(), [](const gene& lhs, const gene& rhs){
        return lhs.fitness < rhs.fitness;
    });
    pool = std::vector<gene>(pool.begin(), pool.begin() + POP_SIZE);
}

void Solver::genetic_algorithm(double time_elapsed_before) {
    parents parent;
    std::vector<gene> offspring;
    auto start = std::chrono::steady_clock::now();
    double time_in_sec;
    double time_left = 570 - time_elapsed_before;
    do {
        for (int j = 0; j < POP_SIZE / 2; ++j) {
            parent_selection(parent);
            crossover(parent, offspring);
            mutation(offspring);
        }
        survivor_selection(offspring);
        offspring.clear();
        std::cout << pool[0].fitness << std::endl;
        auto elapsed = std::chrono::steady_clock::now() - start;
        time_in_sec = std::chrono::duration<double>(elapsed).count();
    } while (time_in_sec < time_left);
    std::cout << "Final HPWL: " << pool[0].fitness << std::endl;
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

void Solver::solve(char *argv[]) {
    auto start = std::chrono::steady_clock::now();
    read(argv);
    create_grid();
    init_pop_from_GP();
    output_file(argv[4]);
    auto elapsed = std::chrono::steady_clock::now() - start;
    double time_in_sec = std::chrono::duration<double>(elapsed).count(); // second in double type
    genetic_algorithm(time_in_sec);
    output_file(argv[4]);
}


