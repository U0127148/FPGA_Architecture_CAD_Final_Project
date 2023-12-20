#include "solver.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_set>

double min_y[4];
double height[4];
int unit[4];
double *region_area;
double **C_matrix_tmp;
double *d_x_tmp, *d_y_tmp;
double **D_matrix, **B_matrix;
double **u, **Z, **D_matrix_inv;

void first(int i, int m, double** arr, double** arr_1) { // 設定前導壹
	double tmp;
	int j;
	tmp = arr[i][i];
	for (j = 0; j < m; j++) {
		arr[i][j] = arr[i][j] / tmp;
		arr_1[i][j] = arr_1[i][j] / tmp;
	}
}

void zero(int i, int m, double** arr, double** arr_1) { // 將前導壹的上下列變為零
	int j, k;
	double tmp;
	for (j = 0; j < m; j++) {
		if (j == i)
			continue;
		tmp = -1 * arr[j][i];
		for (k = 0; k < m; k++) {
			arr[j][k] = arr[j][k] + (tmp * arr[i][k]);
			arr_1[j][k] = arr_1[j][k] + (tmp * arr_1[i][k]);
		}
	}
}

void gauss_jordan(int m, double** arr, double** arr_1) { // 高登-喬登消去法
	int i, j;
	double tmp;
	for (i = 0; i < m; i++) {
		if (arr[i][i] == 0) {
			if (i == (m - 1)) {
				break;
			}
			else {
				for (j = 0; j < m; j++) {
					arr[i][j] = arr[i][j] + arr[i + 1][j];
					arr_1[i][j] = arr_1[i][j] + arr_1[i + 1][j];
				}
				first(i, m, arr, arr_1);
				zero(i, m, arr, arr_1);
			}
		}
		else {
			first(i, m, arr, arr_1);
			zero(i, m, arr, arr_1);
		}
	}
}

double matrixCalculation(size_t row, size_t col, double** m1, double** m2, int q) {
    double val = 0;
    for (int i = 0; i < q; ++i) {
        val += (m1[row][i] * m2[i][col]);
    }
    return val;
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
                    mcell_num++;
                } else if (info == "RAM") {
                    block_type = blockType::RAM;
                    mcell_num++;
                } else if (info == "DSP") {
                    block_type = blockType::DSP;
                    mcell_num++;
                } else if (info == "IO") {
                    block_type = blockType::IO;
                    fcell_num++;
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

                // update cell id (original)
                inst_block->cell_id_org = inst.size() - 1;
                inst_block->cell_id = inst.size() - 1;
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
   
    int cnt = 0;
    for (auto &it : inst) {
        if (it->type == blockType::IO) cnt++;
        else it->cell_id -= cnt;
    }
    // for (auto &it : inst) {
    //     std::cout << it->cell_id_org << ' ' << it->cell_id << std::endl;
    // }
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
    C_matrix = new double*[mcell_num];
    for (size_t i = 0; i < mcell_num; ++i) {
        C_matrix[i] = new double[mcell_num];
    }
    d_x = new double[mcell_num];
    d_y = new double[mcell_num];

    for (size_t i = 0; i < mcell_num; ++i) {
        d_x[i] = 0;
        d_y[i] = 0;
        for (size_t j = 0; j < mcell_num; ++j) {
            C_matrix[i][j] = 0;
        }
    }

    double e;
    Block* block;
    for (auto &net : net_vec) {
        e = 2.0 / (double)net->inst_vec.size();
        int c, c_lamda;
        for (size_t i = 0; i < net->inst_vec.size(); ++i) {
            block = net->inst_vec[i];
            if (block->type == blockType::IO) continue;
            c = block->cell_id;
            C_matrix[c][c] += e * (net->inst_vec.size() - 1);

            for (size_t j = 0; j < net->inst_vec.size(); ++j) {
                if (net->inst_vec[j] == block) continue;
                c_lamda = net->inst_vec[j]->cell_id;
                if (net->inst_vec[j]->type == blockType::IO) {
                    d_x[c] -= e * net->inst_vec[j]->center_x;
                    d_y[c] -= e * net->inst_vec[j]->center_y;
                } else {
                    C_matrix[c][c_lamda] -= e;
                }
            }
        }
    }

    for (auto &cell : inst) {
        cell->size = height[cell->type];
    }

    // reserve space for D matrix and B matrix
    D_matrix = new double*[mcell_num];
    for (size_t i = 0; i < mcell_num; ++i) {
        D_matrix[i] = new double[mcell_num];
    }
    B_matrix = new double*[mcell_num];
    for (size_t i = 0; i < mcell_num; ++i) {
        B_matrix[i] = new double[mcell_num];
    }
    region_area = new double[mcell_num];

    C_matrix_tmp = new double*[mcell_num];
    for (size_t i = 0; i < mcell_num; ++i) {
        C_matrix_tmp[i] = new double[mcell_num];
    }
    d_x_tmp = new double[mcell_num];
    d_y_tmp = new double[mcell_num];

    D_matrix_inv = new double*[mcell_num];
    for (size_t i = 0; i < mcell_num; ++i) {
        D_matrix_inv[i] = new double[mcell_num];
    }
    Z = new double*[mcell_num];
    for (size_t i = 0; i < mcell_num; ++i) {
        Z[i] = new double[mcell_num];
    }
}

void Solver::setUpConstraint() {
    memset(region_area, 0, sizeof(double) * region_num);
    for (auto &cell : inst) {
        region_area[cell->region_idx] += cell->size;
    }

    std::vector<std::vector<Block*>> region_info(region_num);
    for (auto &cell : inst) {
        region_info[cell->region_idx].push_back(cell);
    }

    // values for D matrix
    for (size_t i = 0; i < region_num; ++i) { // row
        for (size_t j = 0; j < region_num; ++j) { // col
            D_matrix[i][j] = 0;
        }
        // the first cell of region_info[i] -> form a diagonal matrix
        D_matrix[i][i] = region_info[i].front()->size / region_area[i];
    }
    int idx = 0;
    // values for B matrix
    for (size_t i = 0; i < region_num; ++i) {
        for (size_t j = 0; j < (mcell_num - region_num); ++j) {
            B_matrix[i][j] = 0;
        }
        // remaining cells of region_info[i]
        for (size_t j = 1; j < region_info[i].size(); ++j) {
            B_matrix[i][idx++] = region_info[i][j]->size / region_area[i];
        }
    }
    // std::cout << "D_matrix: " << std::endl;
    // for (size_t i = 0; i < region_num; ++i) {
    //     for (size_t j = 0; j < region_num; ++j) {
    //         std::cout << D_matrix[i][j] << ' ';
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "B_matrix: " << std::endl;
    // for (size_t i = 0; i < region_num; ++i) {
    //     for (size_t j = 0; j < (mcell_num - region_num); ++j) {
    //         std::cout << B_matrix[i][j] << ' ';
    //     }
    //     std::cout << std::endl;
    // }

    // values for d_x_tmp d_y_tmp C_matrix_tmp
    std::vector<Block*> v;
    for (size_t i = 0; i < region_num; ++i) {
        v.push_back(region_info[i].front());
    }
    for (size_t i = 0; i < region_num; ++i) {
        for (size_t j = 1; j < region_info[i].size(); ++j) {
           v.push_back(region_info[i][j]);
        }
    }
    for (size_t i = 0; i < mcell_num; ++i) {
        d_x_tmp[i] = d_x[v[i]->cell_id];
        d_y_tmp[i] = d_y[v[i]->cell_id];
    }

    for (size_t i = 0; i < mcell_num; ++i) {
        for (size_t j = 0; j < mcell_num; ++j) {
            C_matrix_tmp[i][j] = 0;
            D_matrix_inv[i][j] = 0;
        }
        D_matrix_inv[i][i] = 1;
    }
    for (size_t i = 0; i < mcell_num; ++i) {
        for (size_t j = 0; j < mcell_num; ++j) {
            C_matrix_tmp[i][j] = C_matrix[v[i]->cell_id][v[j]->cell_id];
        }
    }
    gauss_jordan(region_num, D_matrix, D_matrix_inv);

    // calculate Z matrix
    // first calculate -invD*B (region_num, mcell_num-region_num)
    for (size_t i = 0; i < mcell_num; ++i) {
        for (size_t j = 0; j < (mcell_num - region_num); ++j) {
            Z[i][j] = 0;
        }
    }
    for (size_t i = 0; i < region_num; ++i) {
        for (size_t j = 0; j < (mcell_num - region_num); ++j) {
            Z[i][j] = -matrixCalculation(i, j, D_matrix_inv, B_matrix, region_num);
        }
    }
    for (size_t i = region_num; i < (mcell_num - region_num); ++i) {
        Z[i][i] = 1;
    }
}

void Solver::globalOptimize() {

}

void Solver::iterPlacePartition() {
    // while (1) { // need to change k constraint to time constraint

    // }
}