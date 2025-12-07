#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP


#include <vector>

struct Vertex {
    Vertex (): x(0.0f), y(0.0f), z(0.0f), 
               potential(0.0f), density(0.0f), v_id(0) {}

    Vertex (float x, float y, float z): 
           x(x), y(y), z(z), potential(0.0f), density(0.0f), v_id(0) {}

    Vertex (float x, float y, float z, int v_id): 
           x(x), y(y), z(z), potential(0.0f), density(0.0f), v_id(v_id) {}

    Vertex (float x, float y, float z, float potential, int v_id): 
           x(x), y(y), z(z), potential(potential), density(0.0f), v_id(v_id) {}

    Vertex (float x, float y, float z, float potential, float density, int v_id): 
           x(x), y(y), z(z), potential(potential), density(density), v_id(v_id) {}

    Vertex(const Vertex& other): 
           x(other.x), y(other.y), z(other.z),
           potential(other.potential), density(other.density), v_id(other.v_id) {}
    
    Vertex& operator=(const Vertex&) = default;

    bool operator<(const Vertex& other) const {
        return v_id < other.v_id;
    }

    bool operator==(const Vertex& other) const {
        return v_id == other.v_id;
    }    

    float x, y, z;
    float potential;
    float density;
    int v_id;
};

struct Element {
    
    Element()
        : v1(0), v2(0), v3(0), v4(0) {}

    Element(int v1, int v2, int v3, int v4)
        : v1(v1), v2(v2), v3(v3), v4(v4) {}

    Element(int v1, int v2, int v3)
        : v1(v1), v2(v2), v3(v3), v4(0) {}

    Element(const Element& other)
        : v1(other.v1), v2(other.v2), v3(other.v3), v4(other.v4) {}

    Element& operator=(const Element&) = default;

    bool operator==(const Element& other) const {
        return v1 == other.v1 && v2 == other.v2 && v3 == other.v3 && v4 == other.v4;
    }

    int v1, v2, v3, v4; 
};

struct Model {


    int n_vertices;
    int n_elements;
    std::vector<Vertex> vertices;
    std::vector<Element> elements;

    Model() : n_vertices(0), n_elements(0) {}


    Model(const std::vector<Vertex>& verts, const std::vector<Element>& elems)
    : n_vertices(static_cast<int>(verts.size())),
      n_elements(static_cast<int>(elems.size())),
      vertices(verts),
      elements(elems)
    {}



    Model(int n_vertices, int n_elements)
        : n_vertices(n_vertices),
          n_elements(n_elements) {
        vertices.resize(n_vertices);
        elements.resize(n_elements);
    }
    
    void update_counts() {
        n_vertices = static_cast<int>(vertices.size());
        n_elements = static_cast<int>(elements.size());
    }
};



template<typename T_index, typename T_value>
struct CSR_matrix {
    int n_rows;
    int n_cols;
    int n_nonzero;
    std::vector<T_index> row_ptr;
    std::vector<T_index> col_ind;
    std::vector<T_value> values;

    CSR_matrix() : n_rows(0), n_cols(0), n_nonzero(0) {}

    CSR_matrix(int rows, int cols, int nonzero)
        : n_rows(rows), n_cols(cols), n_nonzero(nonzero) {
        row_ptr.resize(n_rows + 1);
        col_ind.resize(n_nonzero);
        values.resize(n_nonzero);
    }
};

#endif