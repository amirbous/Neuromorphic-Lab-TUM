#include <vector>
#include <string>
#include <fstream>

#include "include/model.hpp"

#include <cmath>

void assignBoundaryConditions(Model &model);

template<typename T_index, typename T_value>
void initialize_CSR_indices(const Model &model,
                              CSR_matrix<T_index, T_value> &A);

template<typename T_index, typename T_value>
void fill_FEM_CSR(const Model &model,
                                 const CSR_matrix<T_index, T_value> &A,
                                 std::vector<T_value> &b);


