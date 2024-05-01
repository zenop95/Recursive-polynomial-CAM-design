#include <stdexcept>
#include <fstream>
#include <iostream>
#include <unistd.h>

#include <dace/dace.h>
#include <json/json.h>
#include "cfg/Config.hpp"
#include "astro/AstroRoutines.h"
#include "dynorb/AIDA.h"
#include "geco/JsonParser.hpp"
#include "odeprop/RK78Propagator.hpp"

int main(int argc, char *const argv[])
{
    // Check number of arguments provided
    if (argc < 3)
        throw std::runtime_error("No input or output file specified.");

    // Open desired input file
    std::ifstream inputFile(argv[1], std::ifstream::binary);

    // Read the Json value
    Json::Value root;
    inputFile >> root;

    // Define Json parser
    geco::CJsonParser parser;

    // Load spice kernels
    cfg::EnvConfig().LoadSpiceKernels();

    // Get initial state
    const Json::Value initState = root["initialState"];
    DACE::AlgebraicVector<double> x0(6);
    for (Json::Value::ArrayIndex i = 0; i != initState["state"].size(); i++)
    {
        x0.at(i) = initState["state"][i].asDouble();
    }
    // Read initial epoch
    double et0 = astro::str2et(parser.Read<std::string>(initState, "epoch"));
    // Read mass
    double mass = parser.Read<double>(initState, "mass");
    // Read final epoch
    double etf = astro::str2et(parser.Read<std::string>(root, "finalEpoch"));

    // Get dynamics
    const Json::Value dynamics = root["dynamics"];
    dynorb::AIDA<double> aida(dynamics, parser);

    double A_drag = dynamics["drag"]["area"].asDouble();
    double Cd = dynamics["drag"]["Cd"].asDouble();
    double Bfactor = Cd * A_drag / mass;

    double A_srp = dynamics["srp"]["area"].asDouble();
    double Cr = dynamics["srp"]["Cr"].asDouble();
    double SRPC = Cr * A_srp / mass;

    //Set uncertain parameters
    aida.setUncParam("SRPC", SRPC);
    aida.setUncParam("Bfactor", Bfactor);

    // Perform propagation
    const Json::Value propagJson = root["propagator"];
    odeprop::RK78Propagator<double> rk78(propagJson, parser);
    DACE::AlgebraicVector<double> xf = rk78.propagate(aida, x0, et0, etf);

    // Output state and final epoch
    Json::Value output;
    Json::Value jsonXf(Json::arrayValue);
    for(DACE::AlgebraicVector<double>::const_iterator it=xf.begin(); it!=xf.end(); it++)
    {
        jsonXf.append(Json::Value(*it));
    }
    output["finalState"]["epoch"] = astro::et2str(etf, "YYYY-MM-DDTHR:MN:SC.###### ::UTC");
    output["finalState"]["state"] = jsonXf;
    output["finalState"]["mass"]  = mass;

    std::ofstream outputFile(argv[2]);
    outputFile << output.toStyledString() << std::endl;
}
