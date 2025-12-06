#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP


//struct with only vertices and faces numbers for now
struct ProblemProperties {

	ProblemProperties(): 
	    n_vertices(0), n_faces(0) {}
		
	ProblemProperties(int n_vertices, int n_faces):
	    n_vertices(n_vertices), 
	    n_faces(n_faces) {}

	int n_vertices;
	int n_faces;

};


//struct without load type, just with potential and density --> potential instead of load_value and density
//remove load_type
//we use templates for type of potential and load
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

    float x;
    float y;
    float z;

    float potential;
    float density;
    
    int v_id;
};


struct Element {
    
    Element()
        : v1(0), v2(0), v3(0), v4(0)
          {}


    Element(int v1, int v2, int v3, int v4)
        : v1(v1), v2(v2), v3(v3),v4(v4) 
          {}

    Element(int v1, int v2, int v3)
        : v1(v1), v2(v2), v3(v3), v4(v4)
          {}


    Element(const Element& other)
        : v1(other.v1), v2(other.v2), v3(other.v3), v4(other.v4)
          {}

    Element& operator=(const Element&) = default;

    bool operator==(const Element& other) const {
        return v1 == other.v1 && v2 == other.v2 && v3 == other.v3 && v4 == other.v4;
    }

    int v1;
    int v2;
    int v3;
    int v4;
    
};

#endif