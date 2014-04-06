#pragma once
#include <types.h>
#include <yaml-cpp/yaml.h>
#include <iostream>
#include <string>

namespace YAML
{
	template<>
	struct convert<Grid>{
		static Node encode(const Grid& g){
			Node node;
			node["N_x"] = g.N_x;
			node["N_y"] = g.N_y;
			node["L_x"] = g.L_x;
			node["L_y"] = g.L_y;

			return node;
		}

		static bool decode(const Node& node, Grid& g){
			g.N_x = node["N_x"].as<int>();
			g.N_y = node["N_y"].as<int>();
			g.L_x = node["L_x"].as<Real>();
			g.L_y = node["L_y"].as<Real>();
			return true;
		}
	};

	template<>
	struct convert<Params>{
		static Node encode(const Params& p){
			Node node;
			node["nsteps"] = p.nsteps;
			node["dt"] = p.dt;
			node["Re"] = p.Re;

			return node;
		}

		static bool decode(const Node& node, Params& p){
			p.nsteps = node["nsteps"].as<int>();
			p.dt = node["dt"].as<Real>();
			p.Re = node["Re"].as<Real>();

			return true;
		}
	};



	template<>
	struct convert<Boundaries>{

		static Node encode(const Boundaries& b){
			Node node;

			if (b.left.bc_type == DIRICHLET)
				node["left"][0] = "dirichlet";
			else if (b.left.bc_type == NEUMANN)
				node["left"][0] = "neumann";

			node["left"][1] = b.left.bc_value;


			if (b.right.bc_type == DIRICHLET)
				node["right"][0] = "dirichlet";
			else if (b.right.bc_type == NEUMANN)
				node["right"][0] = "neumann";

			node["right"][1] = b.right.bc_value;


			if (b.bottom.bc_type == DIRICHLET)
				node["bottom"][0] = "dirichlet";
			else if (b.bottom.bc_type == NEUMANN)
				node["bottom"][0] = "neumann";

			node["bottom"][1] = b.left.bc_value;


			if (b.top.bc_type == DIRICHLET)
				node["top"][0] = "dirichlet";
			else if (b.top.bc_type == NEUMANN)
				node["top"][0] = "neumann";

			node["top"][1] = b.left.bc_value;

			return node;
		}

		static bool decode(const Node& node, Boundaries& b){
			// left
			
			if (node["left"][0].as<std::string>() == "dirichlet")
				b.left.bc_type = DIRICHLET;
			else if (node["left"][0].as<std::string>() == "neumann")
				b.left.bc_type = NEUMANN;

			b.left.bc_value = node["left"][1].as<Real>();

			// right
			if (node["right"][0].as<std::string>() == "dirichlet")
				b.right.bc_type = DIRICHLET;
			else if (node["right"][0].as<std::string>() == "neumann")
				b.right.bc_type = NEUMANN;

			b.right.bc_value = node["right"][1].as<Real>();

			// bottom
			if (node["bottom"][0].as<std::string>() == "dirichlet")
				b.bottom.bc_type = DIRICHLET;
			else if (node["bottom"][0].as<std::string>() == "neumann")
				b.bottom.bc_type = NEUMANN;

			b.bottom.bc_value = node["bottom"][1].as<Real>();

			// top
			if (node["top"][0].as<std::string>() == "dirichlet")
				b.top.bc_type = DIRICHLET;
			else if (node["top"][0].as<std::string>() == "neumann")
				b.top.bc_type = NEUMANN;

			b.top.bc_value = node["top"][1].as<Real>();

			return true;

		}
	};
}	// end YAML

