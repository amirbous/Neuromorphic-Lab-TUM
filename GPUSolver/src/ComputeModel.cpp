#include <vector>
#include <set>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <array>
#include <algorithm>
#include <map>

#include "include/model.hpp"
#include "include/ComputeModel.hpp"


/*************************************

Problem: 3d Poisson:  -∇²φ = f  in Ω
                         φ = g  on ∂Ω

            choose f = 40
            φ = 5 * (1 - x² - y² + z²) on ∂Ω (sphere of radius 1)
            f = 30
            g = 5 * (1 - x² - y² + z²) (source orignating from center in a spherical domain)

*******************************************/


/**************************************************
 * FEM first order tetrahedral elements derivation
 * N = 
 * 
 * *************************************************/
 

template<typename T_value>
T_value analytical_solution(T_value x, T_value y, T_value z) {
    return 5.0f * (1.0f - x * x - y * y - z * z);
}


template<typename T_index, typename T_value>
std::vector<T_index> extract_boundary_nodes(const Model<T_index, T_value>& model) {
    
    std::map<Face<T_index>, T_index> face_counts;

    for (const auto& elem : model.elements) {


        T_index element_v[4] = {elem.v1, elem.v2, elem.v3, elem.v4};
        

        //For reference all possible face combinations in a tetrahedron (indexing starts at 1)
        //    {1, 2, 3}, {1, 2, 4}, {1, 3, 4}, {2, 3, 4}



        Face<T_index> f1 = Face<T_index>(elem.v1, elem.v2, elem.v3);
        Face<T_index> f2 = Face<T_index>(elem.v1, elem.v2, elem.v4);
        Face<T_index> f3 = Face<T_index>(elem.v1, elem.v3, elem.v4);
        Face<T_index> f4 = Face<T_index>(elem.v2, elem.v3, elem.v4);

        f1.sort_faces();
        f2.sort_faces();
        f3.sort_faces();
        f4.sort_faces();

        face_counts[f1]++;
        face_counts[f2]++;
        face_counts[f3]++;
        face_counts[f4]++;

    }

    std::set<T_index> boundary_set;
    for (const auto& face_frequency_pair : face_counts) {
        if (face_frequency_pair.second == 1) {
            boundary_set.insert(face_frequency_pair.first.v1);
            boundary_set.insert(face_frequency_pair.first.v2);
            boundary_set.insert(face_frequency_pair.first.v3);
        }
    }

    return std::vector<T_index>(boundary_set.begin(), boundary_set.end());
}

template<typename T_index, typename T_value>
void initialize_CSR_indices(const Model<T_index, T_value> &model,CSR_matrix<T_index, T_value> &A){

    T_index N = model.n_vertices;
    

    A.n_rows = N;
    A.n_cols = N; 
    
    std::vector<std::set<T_index>> adjacency(N);
    // parallelization: ensure if row_node is unique ??
    for (const auto& elem : model.elements) {
        std::vector<T_index> nodes = {elem.v1, elem.v2, elem.v3, elem.v4};
        for (T_index row_node : nodes) { 
            for (T_index col_node : nodes) {

                adjacency[row_node].insert(col_node);
            }
        }
    }

    // room for parallelization openmp : we will see
    T_index total_nnz = 0;
    for (const auto& row_set : adjacency) {
        total_nnz += row_set.size();
    }

    A.n_nonzero = total_nnz;
    A.row_ptr.resize(N + 1);
    A.col_ind.resize(total_nnz);
    A.values.resize(total_nnz, 0.0f); 

    T_index current_nnz_index = 0;
    A.row_ptr[0] = 0;

    for (T_index i = 0; i < N; ++i) {
        
        for (T_index col : adjacency[i]) {
            A.col_ind[current_nnz_index] = col;
            current_nnz_index++;
        }

        A.row_ptr[i + 1] = current_nnz_index;
    }
}





template<typename T_index, typename T_value>
void fill_FEM_CSR(const Model<T_index, T_value> &model,
                                 const CSR_matrix<T_index, T_value> &A,
                                 std::vector<T_value> &b) {

    T_index N = model.n_vertices;
    b.resize(N, 0.0f);

    // will be used later to exclude the boundary nodes from the system
    std::vector<T_index> boundary_nodes = extract_boundary_nodes<T_index, T_value>(model);
    std::cout << "Number of boundary nodes: " << boundary_nodes.size() << std::endl;
    //room for parallelization openmp : we will see
    for (const auto& elem : model.elements) {
        std::vector<int> nodes = {elem.v1, elem.v2, elem.v3, elem.v4};

        T_value x1 = model.vertices[elem.v1].x;
        T_value x2 = model.vertices[elem.v2].x;
        T_value x3 = model.vertices[elem.v3].x;
        T_value x4 = model.vertices[elem.v4].x;

        T_value y1 = model.vertices[elem.v1].y;
        T_value y2 = model.vertices[elem.v2].y;
        T_value y3 = model.vertices[elem.v3].y;
        T_value y4 = model.vertices[elem.v4].y;

        T_value z1 = model.vertices[elem.v1].z;
        T_value z2 = model.vertices[elem.v2].z;
        T_value z3 = model.vertices[elem.v3].z;
        T_value z4 = model.vertices[elem.v4].z;

        T_value volume = fabs(
            (x2 - x1) * ((y3 - y1) * (z4 - z1) - (y4 - y1) * (z3 - z1)) -
            (x3 - x1) * ((y2 - y1) * (z4 - z1) - (y4 - y1) * (z2 - z1)) +
            (x4 - x1) * ((y2 - y1) * (z3 - z1) - (y3 - y1) * (z2 - z1))
        ) / 6.0f;
        // assembling the finite element: this hurts my eyes just to read!!!!!!!!!
        // TODO: actually write the model T_T

    }
}

// Explicit template instantiations to ensure code is emitted for the types we use.

template void initialize_CSR_indices<int, float>(const Model<int, float> &model, CSR_matrix<int, float> &A);
template void initialize_CSR_indices<int, double>(const Model<int, double> &model, CSR_matrix<int, double> &A);

template void fill_FEM_CSR<int, float>(const Model<int, float> &model, const CSR_matrix<int, float> &A,
                                       std::vector<float> &b);
template void fill_FEM_CSR<int, double>(const Model<int, double> &model, const CSR_matrix<int, double> &A,
                                        std::vector<double> &b);

template float analytical_solution<float>(float x, float y, float z);
template double analytical_solution<double>(double x, double y, double z);

template std::vector<int> extract_boundary_nodes<int, float>(const Model<int, float>& model);
template std::vector<int> extract_boundary_nodes<int, double>(const Model<int, double>& model);