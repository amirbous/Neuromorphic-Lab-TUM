
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

            φ = 5 * (1 - x² - y² - z²) on ∂Ω
            f = 30
            g = 5 * (1 - x² - y² - z²)

            ∂Ω is the outside surface of complex volumetric mesh
            Ω is the inner volume of the mesh

*******************************************/


/**************************************************
 * FEM first order tetrahedral elements derivation
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

template<typename T_index>
T_index get_csr_index(const std::vector<T_index>& row_ptr, 
                      const std::vector<T_index>& col_ind, 
                      T_index row, T_index col) {
    
    T_index start = row_ptr[row];
    T_index end = row_ptr[row + 1];
    
    // col_ind is sorted because it was built using std::set iteration
    auto it = std::lower_bound(col_ind.begin() + start, col_ind.begin() + end, col);
    
    if (it != col_ind.begin() + end && *it == col) {
        return (T_index)(it - col_ind.begin());
    }
    return -1; // Should not happen if sparsity pattern is correct
}

template<typename T_index, typename T_value>
void fill_FEM_CSR(const Model<T_index, T_value> &model,
                  CSR_matrix<T_index, T_value> &A,
                  std::vector<T_value> &b) {

    T_index N = model.n_vertices;
    b.assign(N, 0.0f); 

    // 1. ASSEMBLE STIFFNESS MATRIX AND LOAD VECTOR
    for (const auto& elem : model.elements) {
        
        T_index nodes[4] = {elem.v1, elem.v2, elem.v3, elem.v4};
        
        // Cache coordinates per element
        T_value x[4], y[4], z[4];
        for(int i=0; i<4; ++i) {
            x[i] = model.vertices[nodes[i]].x;
            y[i] = model.vertices[nodes[i]].y;
            z[i] = model.vertices[nodes[i]].z;
        }

        // Grad N_i = (Face_Vector_Area) / (3 * Volume)
        
        
        // Standard FEM coeff formulas for Tetrahedrons:
        // b1 = y2(z3-z4) + y3(z4-z2) + y4(z2-z3)
        // c1 = x2(z4-z3) + x3(z2-z4) + x4(z3-z2)
        // d1 = x2(y3-y4) + x3(y4-y2) + x4(y2-y3)
        
        // N0 uses 1,2,3
        // N1 uses 2,0,3
        // V = det([1 x y z]) / 6
        
        T_value B[4], C[4], D[4];

        
        // Re-do robustly using differences to avoid sign confusion
        T_value x21 = x[1]-x[0], x31 = x[2]-x[0], x41 = x[3]-x[0];
        T_value y21 = y[1]-y[0], y31 = y[2]-y[0], y41 = y[3]-y[0];
        T_value z21 = z[1]-z[0], z31 = z[2]-z[0], z41 = z[3]-z[0];
        
        // determinant of tetrahedron matrix
        T_value detJ = x21 * (y31 * z41 - y41 * z31) -
                       x31 * (y21 * z41 - y41 * z21) +
                       x41 * (y21 * z31 - y31 * z21);
                       
        T_value vol = std::abs(detJ) / 6.0f;
        
        //detect degenrate elements
        if (vol < 1e-12) continue;

        // Gradients calculation (dN/dx = b/detJ, dN/dy = c/detJ, dN/dz = d/detJ)
        // We calculate b,c,d directly from definition of adjugate matrix rows
        
        // Node 0 (face 1-3-2) -> indices 1,2,3
        B[0] =   y[1]*(z[2]-z[3]) + y[2]*(z[3]-z[1]) + y[3]*(z[1]-z[2]);
        C[0] = -(x[1]*(z[2]-z[3]) + x[2]*(z[3]-z[1]) + x[3]*(z[1]-z[2]));
        D[0] =   x[1]*(y[2]-y[3]) + x[2]*(y[3]-y[1]) + x[3]*(y[1]-y[2]);
        
        // Node 1 (face 0-2-3) -> indices 0,2,3 (Order 0,3,2 for orientation) 
        // b_i = - det(matrix without row i).

        B[1] =   y[0]*(z[3]-z[2]) + y[2]*(z[0]-z[3]) + y[3]*(z[2]-z[0]);
        C[1] = -(x[0]*(z[3]-z[2]) + x[2]*(z[0]-z[3]) + x[3]*(z[2]-z[0]));
        D[1] =   x[0]*(y[3]-y[2]) + x[2]*(y[0]-y[3]) + x[3]*(y[2]-y[0]);

        // Node 2 (indices 0,1,3)
        B[2] =   y[0]*(z[1]-z[3]) + y[1]*(z[3]-z[0]) + y[3]*(z[0]-z[1]);
        C[2] = -(x[0]*(z[1]-z[3]) + x[1]*(z[3]-z[0]) + x[3]*(z[0]-z[1]));
        D[2] =   x[0]*(y[1]-y[3]) + x[1]*(y[3]-y[0]) + x[3]*(y[0]-y[1]);

        // Node 3 (indices 0,1,2)
        B[3] =   y[0]*(z[2]-z[1]) + y[1]*(z[0]-z[2]) + y[2]*(z[1]-z[0]);
        C[3] = -(x[0]*(z[2]-z[1]) + x[1]*(z[0]-z[2]) + x[2]*(z[1]-z[0]));
        D[3] =   x[0]*(y[2]-y[1]) + x[1]*(y[0]-y[2]) + x[2]*(y[1]-y[0]);

        // Element Load (f = 30)
        T_value f_e = 30.0f;
        T_value load_val = f_e * vol / 4.0f;

        // Assembly
        for (int i = 0; i < 4; ++i) {
            // Add load
            b[nodes[i]] += load_val;

            for (int j = 0; j < 4; ++j) {

                // local stiffness matrix coefficients 
                // K_ij = (grad N_i . grad N_j) * Volume
                // grad N_i = [Bi, Ci, Di] / detJ
                // K_ij = (Bi*Bj + Ci*Cj + Di*Dj) / (detJ * detJ) * vol
                // Since vol = abs(detJ)/6, simplify:
                // K_ij = (Bi*Bj + Ci*Cj + Di*Dj) / (6 * detJ * sign(detJ)) 
                //      = (Bi*Bj + Ci*Cj + Di*Dj) / (6 * abs(detJ))
                
                T_value k_val = (B[i]*B[j] + C[i]*C[j] + D[i]*D[j]) / (6.0f * std::abs(detJ));

                T_index idx = get_csr_index(A.row_ptr, A.col_ind, nodes[i], nodes[j]);
                if (idx != -1) {
                    A.values[idx] += k_val;
                }
            }
        }
    }

    // 2. APPLY BOUNDARY CONDITIONS - rows of boundadry points to unit vectors 
    
    std::vector<T_index> boundary_nodes = extract_boundary_nodes<T_index, T_value>(model);
    
    for (T_index b_node : boundary_nodes) {
        // Calculate g boundary value
        T_value bx = model.vertices[b_node].x;
        T_value by = model.vertices[b_node].y;
        T_value bz = model.vertices[b_node].z;
        T_value g_val = analytical_solution(bx, by, bz);

        // Modify RHS
        b[b_node] = g_val;

        // Modify Matrix A: Row b_node -> [0 0 ... 1 ... 0]
        T_index start = A.row_ptr[b_node];
        T_index end = A.row_ptr[b_node + 1];

        for (T_index k = start; k < end; ++k) {
            if (A.col_ind[k] == b_node) {
                A.values[k] = 1.0f; // Diagonal
            } else {
                A.values[k] = 0.0f; // Off-diagonal
            }
        }
    }
}


// Explicit template instantiations to ensure code is emitted for the types we use.

template void initialize_CSR_indices<int, float>(const Model<int, float> &model, CSR_matrix<int, float> &A);
template void initialize_CSR_indices<int, double>(const Model<int, double> &model, CSR_matrix<int, double> &A);

template void fill_FEM_CSR<int, float>(const Model<int, float> &model,
                                      CSR_matrix<int, float> &A,
                                      std::vector<float> &b);
template void fill_FEM_CSR<int, double>(const Model<int, double> &model,
                                       CSR_matrix<int, double> &A,
                                       std::vector<double> &b);

template float analytical_solution<float>(float x, float y, float z);
template double analytical_solution<double>(double x, double y, double z);

template std::vector<int> extract_boundary_nodes<int, float>(const Model<int, float>& model);
template std::vector<int> extract_boundary_nodes<int, double>(const Model<int, double>& model);

template int get_csr_index(const std::vector<int>& row_ptr, 
                      const std::vector<int>& col_ind, 
                      int row, int col);
