#include "solver.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_set>

double min_y[4];
double height[4];
int unit[4];
std::unordered_map<Net*, std::vector<Net*>> adj;

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
                Block* resource_block = new Block(block_type, center_x, center_y);
                resource.push_back(resource_block);
                nameToResource[name] = resource_block;
            }
            n++;
        }
    }
    fin.close();

    // second input file (inst)
    fin.open(argv[2]);
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
                if (info == "CLB") {
                    block_type = blockType::CLB;
                } else if (info == "RAM") {
                    block_type = blockType::RAM;
                } else if (info == "DSP") {
                    block_type = blockType::DSP;
                } else if (info == "IO") {
                    block_type = blockType::IO;
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
                Block* inst_block = new Block(block_type, center_x, center_y);
                inst.push_back(inst_block);
                nameToInst[name] = inst_block;
            }
            n++;
        }
    }
    fin.close();

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
                net->inst_vec.push_back(nameToInst[info]);
            }
            n++;
        }
        net_vec.push_back(net);
    }
    fin.close();

    /*
    for (auto &net : net_vec) {
        net->init();
        // std::cout << net->name << " HPWL: " << net->HPWL << std::endl;
    }
    */
}

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
    
    /*
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
    */

    std::vector<Block*> single_window[4];
    int u = resource[0]->unit;
    for (auto &it : resource) {
        if (it->unit != u) {
            window[blockType::CLB].push_back(single_window[blockType::CLB]);
            window[blockType::RAM].push_back(single_window[blockType::RAM]);
            window[blockType::DSP].push_back(single_window[blockType::DSP]);
            single_window[blockType::CLB].clear();
            single_window[blockType::RAM].clear();
            single_window[blockType::DSP].clear();
            u = it->unit;
            single_window[it->type].push_back(it);
        } else {
            single_window[it->type].push_back(it);
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

void Solver::setUpObject() {

}