#include <vector>
#include <set>
#include <string>
#include <fstream>

#include "math.h"

#include "include/model.hpp"
#include "include/ComputeModel.hpp"

void assignBoundaryConditions(Model &model){
    for (int i = 0; i < model.n_vertices; ++i) {
        model.vertices[i].potential = sqrt(model.vertices[i].x * model.vertices[i].x + 
                                           model.vertices[i].y * model.vertices[i].y + 
                                           model.vertices[i].z * model.vertices[i].z);
    }
}


template<typename T_index, typename T_value>
void initialize_CSR_indices(const Model &model,CSR_matrix<T_index, T_value> &A){

    T_index N = model.n_vertices;
    

    A.n_rows = N;
    A.n_cols = N; 
    
    std::vector<std::set<int>> adjacency(N);


    for (const auto& elem : model.elements) {
        std::vector<int> nodes = {elem.v1, elem.v2, elem.v3, elem.v4};
        for (T_index row_node : nodes) { 
            for (T_index col_node : nodes) {

                adjacency[row_node].insert(col_node);
            }
        }
    }

    int total_nnz = 0;
    for (const auto& row_set : adjacency) {
        total_nnz += row_set.size();
    }

    A.n_nonzero = total_nnz;
    A.row_ptr.resize(N + 1);
    A.col_ind.resize(total_nnz);
    A.values.resize(total_nnz, 0.0f); 

    int current_nnz_index = 0;
    A.row_ptr[0] = 0;

    for (T_index i = 0; i < N; ++i) {
        
        for (int col : adjacency[i]) {
            A.col_ind[current_nnz_index] = col;
            current_nnz_index++;
        }

        A.row_ptr[i + 1] = current_nnz_index;
    }
}


template<typename T_index, typename T_value>
void fill_FEM_CSR(const Model &model,
                                 const CSR_matrix<T_index, T_value> &A,
                                 std::vector<T_value> &b) {

    T_index N = model.n_vertices;
    b.resize(N, 0.0f);


    for (const auto& elem : model.elements) {
        std::vector<int> nodes = {elem.v1, elem.v2, elem.v3, elem.v4};

        float x1 = model.vertices[elem.v1].x;
        float x2 = model.vertices[elem.v2].x;
        float x3 = model.vertices[elem.v3].x;
        float x4 = model.vertices[elem.v4].x;

        float y1 = model.vertices[elem.v1].y;
        float y2 = model.vertices[elem.v2].y;
        float y3 = model.vertices[elem.v3].y;
        float y4 = model.vertices[elem.v4].y;

        float z1 = model.vertices[elem.v1].z;
        float z2 = model.vertices[elem.v2].z;
        float z3 = model.vertices[elem.v3].z;
        float z4 = model.vertices[elem.v4].z;

        T_value volume = fabs(
            (x2 - x1) * ((y3 - y1) * (z4 - z1) - (y4 - y1) * (z3 - z1)) -
            (x3 - x1) * ((y2 - y1) * (z4 - z1) - (y4 - y1) * (z2 - z1)) +
            (x4 - x1) * ((y2 - y1) * (z3 - z1) - (y3 - y1) * (z2 - z1))
        ) / 6.0f;

    }
}

// Explicit template instantiations to ensure code is emitted for the types we use.
template void initialize_CSR_indices<int, float>(const Model &model, CSR_matrix<int, float> &A);
template void initialize_CSR_indices<int, double>(const Model &model, CSR_matrix<int, double> &A);

template void fill_FEM_CSR<int, float>(const Model &model, const CSR_matrix<int, float> &A,
                                       std::vector<float> &b);
template void fill_FEM_CSR<int, double>(const Model &model, const CSR_matrix<int, double> &A,
                                        std::vector<double> &b);