#include "include/IO.hpp"
#include "include/model.hpp"

#include <iostream>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include <fstream>
#include <string>
#include <iomanip>

void DisplayWelcomeHeader()
{

        std::cout << "**************************************" << std::endl;
        std::cout << "*          3D FEM POISSON            *" << std::endl;
        std::cout << "**************************************" << std::endl;
        std::cout << std::endl;
}


// TODO: adapt for tetrahedral elements - DONE
void ReadVTK(std::string model_name, Model &model)
{

        std::string vtk_filename = model_name + ".vtk";
        std::string vtk_file = "data/" + vtk_filename;
        std::ifstream vtk_file_stream(vtk_file);


        std::vector<Vertex>& vertices = model.vertices;
        std::vector<Element>& elements = model.elements;

        if (!vtk_file_stream.is_open())
        {
                std::cerr << "Error: Could not open VTK file " << vtk_filename << std::endl;
                return;
        }

        std::string buffer;

        std::vector<unsigned long int> offsets;
        std::vector<unsigned long int> connectivity_raw;
        unsigned long int n_cells_header = 0;

        while (vtk_file_stream >> buffer)
        {

                if (buffer == "POINTS")
                {
                        int n_points;
                        std::string type;
                        vtk_file_stream >> n_points >> type;

                        model.n_vertices = n_points;
                        vertices.resize(n_points);

                        for (int i = 0; i < n_points; ++i)
                        {
                                float x, y, z;
                                vtk_file_stream >> x >> y >> z;
                                vertices[i] = Vertex(x, y, z, i);
                        }
                }

                else if (buffer == "CELLS")
                {
                        int size_dummy;
                        vtk_file_stream >> n_cells_header >> size_dummy;
                        n_cells_header --;
                        model.n_elements = n_cells_header;
                        elements.reserve(n_cells_header);
                }


                else if (buffer == "OFFSETS")
                {
                        std::string type;
                        vtk_file_stream >> type;
                        offsets.resize(n_cells_header + 1);

                        for (int i = 0; i <= n_cells_header; ++i)
                        {
                                vtk_file_stream >> offsets[i];
                        }
                }

                else if (buffer == "CONNECTIVITY")
                {
                        std::string type;
                        vtk_file_stream >> type; 


                        unsigned long int total_indices = offsets.back();
                        connectivity_raw.resize(total_indices);
                        std::cout << "Total connectivity indices to read: " << total_indices << std::endl;
                        for (unsigned long int i = 0; i < total_indices; ++i)
                        {
                                vtk_file_stream >> connectivity_raw[i];
                        }


                        for (int i = 0; i < n_cells_header; ++i)
                        {
                                unsigned long int start = offsets[i];
                                unsigned long int end = offsets[i + 1];
                                unsigned long int count = end - start;

                                if (count == 4)
                                {

                                        int v1 = static_cast<int>(connectivity_raw[start]);
                                        int v2 = static_cast<int>(connectivity_raw[start + 1]);
                                        int v3 = static_cast<int>(connectivity_raw[start + 2]);
                                        int v4 = static_cast<int>(connectivity_raw[start + 3]);

                                        elements.emplace_back(v1, v2, v3, v4);
                                }
                                else
                                {
                                        std::cerr << "Offset did not match a tetrahedra element: Cell " << i << " has " << count << " vertices " << std::endl;
                                        return;
                                }
                        }
                        connectivity_raw.clear();
                        connectivity_raw.shrink_to_fit();

                        offsets.clear();
                        offsets.shrink_to_fit();
                }
        }
        
        vtk_file_stream.close();
}


// TODO: adapt for tetrahedral elements - DONE
void write_vtu(const std::string problem_name,
               Model &model)
{
    std::ofstream fstream;
    std::string fname = problem_name + "_output.vtu";

    std::vector<Vertex>& vertices = model.vertices;
    std::vector<Element>& elements = model.elements;

    int nvertices = model.n_vertices;
    int nelements = model.n_elements; 

    fstream.open(fname);
    if (fstream.is_open())
    {
        fstream << "<?xml version=\"1.0\"?>" << std::endl;
        fstream << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"";
        fstream << " byte_order=\"LittleEndian\" header_type=\"UInt64\">" << std::endl;
        fstream << "<UnstructuredGrid>" << std::endl;
        fstream << "<Piece NumberOfPoints=\"" << nvertices << "\"  ";
        fstream << "NumberOfCells=\"" << nelements << "\">" << std::endl;
        

        fstream << "<Points>" << std::endl;
        fstream << "<DataArray type=\"Float32\" Name=\"Points\"";
        fstream << " NumberOfComponents=\"3\" Format=\"ascii\">" << std::endl;

        for (size_t i = 0; i < nvertices; i++)
        {
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

        for (size_t i = 0; i < nelements; i++)
        {

            fstream << elements[i].v1
                    << " " << elements[i].v2
                    << " " << elements[i].v3
                    << " " << elements[i].v4
                    << std::endl;
        }
        fstream << "</DataArray>" << std::endl;

 
        fstream << "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << std::endl;
        
        long long current_offset = 0;
        for (int i = 0; i < nelements; i++)
        {
            current_offset += 4; 
            fstream << current_offset << " ";
            if ((i + 1) % 20 == 0) fstream << std::endl;
        }
        fstream << std::endl;
        fstream << "</DataArray>" << std::endl;

        fstream << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;

        for (size_t i = 0; i < nelements; i++)
        {
            fstream << "10" << " ";
            if ((i + 1) % 20 == 0) fstream << std::endl;
        }
        fstream << std::endl;
        fstream << "</DataArray>" << std::endl;
        fstream << "</Cells>" << std::endl;


        fstream << "<PointData>" << std::endl;

        fstream << "<DataArray type=\"Float32\" Name=\"Potential\" ";
        fstream << "NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
        for (int i = 0; i < nvertices; i++)
        {
            fstream << vertices[i].potential << std::endl;
        }
        fstream << "</DataArray>" << std::endl;


        fstream << "<DataArray type=\"Float32\" Name=\"density\" ";
        fstream << "NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
        for (int i = 0; i < nvertices; i++)
        {
            fstream << vertices[i].density << std::endl;
        }
        fstream << "</DataArray>" << std::endl;

        fstream << "</PointData>" << std::endl;

        // --- Footer ---
        fstream << "</Piece>" << std::endl;
        fstream << "</UnstructuredGrid>" << std::endl;
        fstream << "</VTKFile>";

        fstream.close();
    }
    else
    {
        std::cerr << "Error in opening result file" << std::endl;
    }
}