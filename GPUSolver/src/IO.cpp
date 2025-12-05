#include "include/IO.hpp"

#include <iostream>

#include <fstream>
#include <string>
#include <iomanip>


void DisplayWelcomeHeader() {

    std::cout << "**************************************" << std::endl;
    std::cout << "*          3D FEM POISSON            *" << std::endl;
    std::cout << "**************************************" << std::endl;
    std::cout << std::endl;
}


void ReadVTK(std::string model_name, std::vector<Vertex>& vertices, std::vector<Face>& faces,
				ProblemProperties& problem_properties) {

    std::string vtk_filename = model_name + ".vtk";

    std::ifstream vtk_file_stream(vtk_filename);
    if (!vtk_file_stream.is_open()) {
        std::cerr << "Error: Could not open VTK file " << vtk_filename << std::endl;
        return;
    }

    return;
}


void ReadObj(std::string model_name, std::vector<Vertex>& vertices, std::vector<Face>& faces,
                ProblemProperties& problem_properties) {

    std::string obj_filename = model_name + ".obj";

    std::ifstream obj_file_stream(obj_filename);
    if (!obj_file_stream.is_open()) {
        std::cerr << "Error: Could not open OBJ file " << obj_filename << std::endl;
        return;
    }

    return;
}
    

void write_vtu(const std::string problem_name, const std::vector<Vertex>& vertices, const std::vector<Face>& faces,
                                                        ProblemProperties& problem_properties) {


        std::ofstream fstream;


        std::string fname = problem_name + ".vtu";

        int nvertices = problem_properties.n_vertices;
        int nfaces    = problem_properties.n_faces;

        fstream.open(fname);
        if (fstream.is_open()) {

                // headers
                fstream << "<?xml version=\"1.0\"?>" << std::endl;
                fstream << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"";
                fstream << " byte_order=\"LittleEndian\" header_type=\"UInt64\">" << std::endl;
                fstream << "<UnstructuredGrid>" << std::endl;
                fstream << "<Piece NumberOfPoints=\"" << problem_properties.n_vertices << "\"  ";
                fstream << "NumberOfCells=\"" << problem_properties.n_faces << "\">" << std::endl;
                fstream << "<Points>" << std::endl;
                fstream << "<DataArray type=\"Float32\" Name=\"Points\"";
                fstream << " NumberOfComponents=\"3\" Format=\"ascii\">" << std::endl;

                // Writing the vertices
                for (size_t i = 0; i < problem_properties.n_vertices; i++) {

                        
                        fstream << std::fixed << std::setprecision(7) << vertices[i].x
                         << " " << std::fixed << std::setprecision(7) << vertices[i].y
                         << " " << std::fixed << std::setprecision(7) << vertices[i].z
                         << std::endl;

                }

                fstream << "</DataArray>" << std::endl;
                fstream << "</Points>" << std::endl;
                fstream << "<Cells>" << std::endl;
                fstream << "<DataArray type=\"Int64\" ";
                fstream << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

                for (size_t i = 0; i < nfaces; i++) {
                        fstream << faces[i].v1
                         << " " << faces[i].v2
                         << " " << faces[i].v3
                         << std::endl;
                }

                fstream << "</DataArray>" << std::endl;

                fstream << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
                for (int i = 3; i < nfaces * 3 + 1; i+= 3) {
                        fstream << i << ((i % 30) != 0 ? " " : " \n");
                }

                if (nfaces  *  3 % 30 != 0) {
                        fstream << std::endl;
                }

                fstream << "</DataArray>" << std::endl;


                fstream << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;

                for (size_t i = 0; i < nfaces; i++) {
                        fstream << "5" << ((i+1) % 10 != 0 ? " " : " \n");
                }

                if (nfaces % 10 != 0) {
                        fstream << std::endl;
                }
                fstream << "</DataArray>" << std::endl;
                fstream << "</Cells>" << std::endl;

                fstream << "<PointData>" << std::endl;
                fstream << "<DataArray type=\"Float32\" Name=\"Potential\" ";
                fstream << "NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;

                for (int i = 0; i < nvertices; i++) {
                        fstream << vertices[i].potential << std::endl;
                }

                fstream << "</DataArray>" << std::endl;

                fstream << "<DataArray type=\"Float32\" Name=\"density\" ";
                fstream << "NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
                for (int i = 0; i < nvertices; i++) {
                        fstream << vertices[i].density << std::endl;
                }

                fstream << "</DataArray>" << std::endl;

                fstream << "</PointData>" << std::endl;

                fstream << "</Piece>" << std::endl;
                fstream << "</UnstructuredGrid>" << std::endl;
                fstream << "</VTKFile>";


                fstream.close();


        } else {

                std::cerr << "Error in opening result file" << std::endl;
        }

}
