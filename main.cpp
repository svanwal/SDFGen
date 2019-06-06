// SDFGen - A simple commandline utility to convert polyhedron shape models of small bodies into signed distance fields.
// This code is based on the GitHub repo by Christopher Batty, which was released under the MIT License
// The original header and MIT License are copied at the bottom of this main file.

#include "makelevelset3.h"
#include "config.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>
#include <string>
#include <ctime>

int main(int argc, char* argv[]) {

    // Checking number of inputs
    if(argc<4) {
        std::cout << "******\n";
        std::cout << "SDFGen - A utility for converting polyhedron shape models of small bodies into grid-based signed distance fields.\n";
        std::cout << "******\n\n";
        std::cout << "Usage: ./bin/SDFGen <filename> <dx> <padding> <format>\n";
        std::cout << "  <filename> specifies a Wavefront OBJ (text) file representing a *triangle* mesh (no quad or poly meshes allowed). File must use the suffix \".obj\".\n";
        std::cout << "  <dx> specifies the length of grid cell in the resulting distance field.\n";
        std::cout << "  <padding> specifies the number of cells worth of padding between the object bounding box and the boundary of the distance field grid. Minimum is 1.\n\n";
        std::cout << "   <format> is an optional argument that specifies the output format, use 0 to save as text file (default), 1 to save as binary file, and 2 to save as both.\n";
        std::cout << "The output filename will match that of the input, with the OBJ suffix replaced with SDF. The output file format is:\n";
        std::cout << "  <ni> <nj> <nk> the integer dimensions of the resulting distance field\n";
        std::cout << "  <origin_x> <origin_y> <origin_z> the 3D position of the grid origin\n";
        std::cout << "  <dx> the grid spacing\n";
        std::cout << "  <d000> <d001> <d002> <...> the signed distance data values, in ascending order of i, then j, then k\n\n";
        exit(-1);
    }

    // Reading inputs
    std::string filename(argv[1]);
    if(filename.size() < 5 || filename.substr(filename.size()-4) != std::string(".obj")) {
      std::cerr << "Error: Expected OBJ file with filename of the form <name>.obj.\n";
      exit(-1);
    }
    // Read cell size
    std::stringstream arg2(argv[2]);
    float dx;
    arg2 >> dx;
    // Read padding
    std::stringstream arg3(argv[3]);
    int padding;
    arg3 >> padding;
    if(padding<1){padding = 1;}
    // Read output format
    int format = 0;
    if(argc>4){
        std::stringstream arg4(argv[4]);
        arg4 >> format;
        if(format>2){
            std::cerr << "Error: Output format must be 0, 1, or 2. You entered the value <" << format << ">\n";
            exit(-1);
        }
    }
    
    // Start with a massive inside out bound box.
    Vec3f min_box(std::numeric_limits<float>::max(),std::numeric_limits<float>::max(),std::numeric_limits<float>::max());
    Vec3f max_box(-std::numeric_limits<float>::max(),-std::numeric_limits<float>::max(),-std::numeric_limits<float>::max());

    // Reading vertices and facets
    std::cout << "Reading .OBJ shape model <" << filename << ">...\n";
    std::ifstream infile(argv[1]);
    if(!infile) {
      std::cerr << "   Failed to open. Is it in the /SDFGen/data/ folder? Terminating.\n";
      exit(-1);
    }

    int ignored_lines = 0;
    std::string line;
    std::vector<Vec3f>  pts;
    std::vector<Vec3ui> tri;
    uint f = 0;
    
    while(!infile.eof()){ //https://en.wikipedia.org/wiki/Wavefront_.obj_file
        std::getline(infile, line);
        if(line.substr(0,1) == std::string("v") && line.substr(0,2) != std::string("vn")){ // Reading a vertex
            std::stringstream data(line);
            char c;
            Vec3f point;
            data >> c >> point[0] >> point[1] >> point[2];
            pts.push_back(point);
            update_minmax(point, min_box, max_box);
        }
        else if(line.substr(0,1) == std::string("f")){ // Reading a facet
            std::stringstream data(line);
            char c;
            std::string f1, f2, f3;
            data >> c >> f1 >> f2 >> f3;
            uint l1 = f1.length();
            uint l2 = f2.length();
            uint l3 = f3.length();
            uint k1, k2, k3;
            if(l1<4 || l2<4 || l3<4){ // No texture/color info in the facet
                k1 = atoi(f1.c_str()) - 1;
                k2 = atoi(f2.c_str()) - 1;
                k3 = atoi(f3.c_str()) - 1;
            }
            else{ // There may be texture/color info that needs to be cut 
                std::string m1 = f1.substr((l1-2)/2,2);
                std::string m2 = f2.substr((l2-2)/2,2);
                std::string m3 = f3.substr((l3-2)/2,2);
                if(m1==std::string("//")){ // There is texture/color info
                    k1 = atoi(f1.substr(0,(l1-2)/2).c_str()) - 1;
                }
                else{ // There is not
                    k1 = atoi(f1.c_str()) - 1;
                }
                if(m2==std::string("//")){ // There is texture/color info
                    k2 = atoi(f2.substr(0,(l2-2)/2).c_str()) - 1;
                }
                else{ // There is not
                    k2 = atoi(f2.c_str()) - 1;
                }
                if(m3==std::string("//")){ /// There is texture/color info
                    k3 = atoi(f3.substr(0,(l3-2)/2).c_str()) - 1;
                }
                else{ // There is not
                    k3 = atoi(f3.c_str()) - 1;
                }
            }
            tri.push_back(Vec3ui(k1,k2,k3));
            f++;
        }
        else{
            ++ignored_lines;
        }
    }
    infile.close();
    std::cout << "   Read in " << pts.size() << " vertices and " << tri.size() << " faces." << std::endl;
    if(tri.size()==2*pts.size()-4){
        std::cout << "   This model satisfies Euler's rule for watertight triangular polyhedra.\n";
    }
    else{
        std::cout << "   WARNING: This model DOES NOT satisfy Euler's rule for watertight triangular polyhedra.\n";
        std::cout << "            Based on the number of vertices, it should have " << 2*pts.size()-4 << " facets instead.\n";
        std::cout << "            SDFGen will proceed, but use the resulting SDF at your own risk...\n";
    }
    
    // Add padding around the box.
    std::cout << "Computing SDF...\n";
    Vec3f unit(1,1,1);
    min_box -= padding*dx*unit;
    max_box += padding*dx*unit;
    Vec3ui sizes = Vec3ui((max_box - min_box)/dx);
    std::cout << "   Origin of the SDF: (" << min_box << ")\n";
    std::cout << "   Maximum extent of the SDF: (" << max_box << ")\n";
    std::cout << "   Dimensions of the SDF: (" << sizes << ") cells\n";
    Array3f d_grid;
    std::clock_t t_start = std::clock();
    make_level_set3(tri, pts, min_box, dx, sizes[0], sizes[1], sizes[2], d_grid); // Compute the SDF
    std::clock_t t_end = std::clock();
    std::cout << "\n   The SDF has been computed!\n";
    std::cout << "   Elapsed time: " << (t_end - t_start)/(float)CLOCKS_PER_SEC << " seconds.\n";

    if(format==0 || format==2){
        // Saving as text file
        std::string outname = filename.substr(0, filename.size()-4) + std::string(".txt");
        std::cout << "Saving the SDF to text file <" << outname << ">...\n";
        std::ofstream fs;
        fs.open(outname.c_str());
        fs.precision(8); // Select precision here
        if(!fs.is_open()){
            std::cerr << "   Error: Unable to open output file! Check the permissions.\n";
            exit(-1);
        }
        fs << d_grid.ni << " " << d_grid.nj << " " << d_grid.nk << "\n";
        fs << min_box[0] << " " << min_box[1] << " " << min_box[2] << "\n";
        fs << dx << "\n";
        uint p = 0;
        for(uint k=0; k<d_grid.nk; k++){
            std::cout << "\r   Progress " << 100.0*(float)k/(float)(d_grid.nk-1) << "%                ";
            for(uint j=0; j<d_grid.nj; j++){
                for(uint i=0; i<d_grid.ni; i++){
                    fs << d_grid.a[p] << " ";
                    p++;
                }
            }
        }
        fs.close();
    }
    if(format==1 || format==2){
        // Saving as binary file
        std::string outname = filename.substr(0, filename.size()-4) + std::string(".sdf");
        std::cout << "Saving the SDF to binary file <" << outname << ">...\n";
        std::ofstream fs(outname.c_str(), std::ios::binary | std::ios::out);
        if(!fs.is_open()){
            std::cerr << "   Error: Unable to open output file! Check the permissions.\n";
            exit(-1);
        }
        fs.write(reinterpret_cast<const char*>(&d_grid.ni),sizeof(d_grid.ni));
        fs.write(reinterpret_cast<const char*>(&d_grid.nj),sizeof(d_grid.nj));
        fs.write(reinterpret_cast<const char*>(&d_grid.nk),sizeof(d_grid.nk));
        fs.write(reinterpret_cast<const char*>(&min_box[0]),sizeof(min_box[0]));
        fs.write(reinterpret_cast<const char*>(&min_box[1]),sizeof(min_box[0]));
        fs.write(reinterpret_cast<const char*>(&min_box[2]),sizeof(min_box[0]));
        fs.write(reinterpret_cast<const char*>(&dx),sizeof(dx));
        for(unsigned int i = 0; i < d_grid.a.size(); ++i) {
            fs.write(reinterpret_cast<const char*>(&d_grid.a[i]),sizeof(d_grid.a[i]));
        }
        fs.close();
    }
    std::cout << "\n   Saving completed.\n";

    std::cout << "SDFGen has finished. See you next time!\n";
    return 0;
    
}

//SDFGen - A simple grid-based signed distance field (level set) generator for triangle meshes.
//Written by Christopher Batty (christopherbatty@yahoo.com, www.cs.columbia.edu/~batty)
//...primarily using code from Robert Bridson's website (www.cs.ubc.ca/~rbridson)
//This code is public domain. Feel free to mess with it, let me know if you like it.
//
//The MIT License (MIT)
//
//Copyright (c) 2015, Christopher Batty
//
//Permission is hereby granted, free of charge, to any person obtaining a copy of
//this software and associated documentation files (the "Software"), to deal in
//the Software without restriction, including without limitation the rights to
//use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
//of the Software, and to permit persons to whom the Software is furnished to do
//so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
//FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
//COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
//IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
//CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.