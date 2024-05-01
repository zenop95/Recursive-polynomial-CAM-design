#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

template <class T>
void find_field(const std::string &fieldname, std::ifstream &file, T &out)
{
    file.seekg(0, file.beg);

    std::string aux;

    // Model name extraction
    do
    {
        getline(file, aux);
    } while (aux.substr(0, fieldname.length()).compare(fieldname) != 0);

    std::stringstream tmp(aux);

    tmp >> aux >> out;
}

void read_coefficients(std::ifstream &file, int max_deg, double *C, double *S)
{

    int n_elem = max_deg + 1;

    file.seekg(0, file.beg);

    std::string fieldname("gfc");

    std::string aux;
    do
    {
        getline(file, aux);
    } while (aux.substr(0, fieldname.length()).compare(fieldname) != 0);

    int i, j;

    // Iterate over elements
    getline(file, aux);
    while ((aux.substr(0, fieldname.length()).compare(fieldname) == 0) && (file.good()))
    {

        std::stringstream tmp(aux);
        double number;

        tmp >> aux >> i >> j;
        tmp >> number;
        C[i * n_elem + j] = number;
        tmp >> number;
        S[i * n_elem + j] = number;
        
        getline(file, aux);
    }

    return;
}

int main(int argc, char *const argv[])
{
    // Check number of arguments provided
    if (argc == 1)
        throw std::runtime_error("No input file specified.");

    // Open desired file
    std::ifstream gravfile(argv[1]);

    //Check if file exists
    if (!gravfile.is_open())
        throw std::runtime_error("File not found.");

    // Create binary file to store values
    std::string binfname(argv[1]);
    size_t pos = binfname.find_last_of(".");
    binfname = binfname.substr(0, pos) + ".bin";
    if (argc > 2)
        binfname = argv[2];

    std::ofstream bingravfile(binfname.c_str(), std::ostream::binary);

    // Extract model name
    std::string name;
    find_field("modelname", gravfile, name);

    // Exctract gravitational parameter
    double GM;
    find_field("earth_gravity_constant", gravfile, GM);
    GM = GM * 1e-9;

    // Extract radius
    double radius;
    find_field("radius", gravfile, radius);
    radius = radius * 1e-3;

    // Extract max degree
    int max_deg;
    find_field("max_degree", gravfile, max_deg);

    double *C = new double[(max_deg + 1) * (max_deg + 1)];
    double *S = new double[(max_deg + 1) * (max_deg + 1)];

    read_coefficients(gravfile, max_deg, C, S);

    // Write extracted values to binary file
    bingravfile.write((char *)&max_deg, sizeof(int));
    bingravfile.write((char *)&GM, sizeof(double));
    bingravfile.write((char *)&radius, sizeof(double));

    bingravfile.write((char *)C, sizeof(double) * (max_deg + 1) * (max_deg + 1));
    bingravfile.write((char *)S, sizeof(double) * (max_deg + 1) * (max_deg + 1));

    gravfile.close();
    bingravfile.close();
    delete[] C;
    delete[] S;
    return 0;
}
