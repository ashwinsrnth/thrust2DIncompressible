#include <read.h>
#include <parser.h>

void read_inputs(std::string file_name, Grid& grid, Params& params, Boundaries& boundaries){
    /*
        Reads the input file specified by file_name 
    */

    YAML::Node node = YAML::LoadFile("simdata.yaml");

    grid = node["grid"].as<Grid>();
    params = node["params"].as<Params>();
    boundaries = node["boundaries"].as<Boundaries>();

}