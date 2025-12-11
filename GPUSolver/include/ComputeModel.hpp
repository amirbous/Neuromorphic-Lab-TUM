#include <vector>
#include <string>
#include <fstream>

#include "include/model.hpp"

#include <cmath>


template<typename T_index, typename T_value>
void initialize_CSR_indices(const Model<T_index, T_value> &model,
                              CSR_matrix<T_index, T_value> &A);

template<typename T_index, typename T_value>
void fill_FEM_CSR(const Model<T_index, T_value> &model,
                                 const CSR_matrix<T_index, T_value> &A,
                                 std::vector<T_value> &b);

template<typename T_index>
std::vector<T_index> extract_boundary_nodes(const Model<T_index, double>& model);

template<typename T_value>
T_value analytical_solution(T_value x, T_value y, T_value z);



template<typename T_index, typename T_value>
void fill_FEM_CSR(const Model<T_index, T_value> &model,
                  CSR_matrix<T_index, T_value> &A,
                  std::vector<T_value> &b);
        

