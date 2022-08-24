#pragma once
#include "Common.h"
#include "Json.h"
#include <fstream>

namespace Meso {
	class MetaData
	{
	public:
		std::string output_base_dir;
		bf::path base_path;
		std::ofstream data_output;
	};

	class DriverMetaData : public MetaData
	{
	public:
		//fixed part
		int fps = 25;
		real cfl = 1.0;
		real time_per_frame = 0.04;
		real min_step_frame_fraction = 0;	//if set to 0.1, it means the minimal iteration time is 0.1*time_per_frame

		int first_frame;
		int last_frame;
		
		//fill for every time step
		int current_frame;
		real current_time;
		real dt;

		real Time_At_Frame(int frame) {
			return frame * time_per_frame;
		}

		void Init(json& j) {
			output_base_dir = Json::Value(j, "output_base_dir", std::string("output"));
			base_path = bf::path(output_base_dir);

			fps = Json::Value(j, "fps", 25);
			cfl = Json::Value(j, "cfl", (real)1.0);
			time_per_frame = 1.0 / fps;
			min_step_frame_fraction = Json::Value(j, "min_step_frame_fraction", (real)0);

			first_frame = Json::Value(j, "first_frame", 0);
			last_frame = Json::Value(j, "last_frame", fps * 10);
			
			current_frame = first_frame;
			current_time = Time_At_Frame(current_frame);
			dt = 1.0 / fps;
		}
	};

	class OptimizerDriverMetaData : public MetaData
	{
	public:
		real tol;
		int iter_count;
		int max_iter_num;

		void Init(json& j) {
			tol = Json::Value(j, "tol", (real) 1e-7);
			max_iter_num = Json::Value(j, "max_iter_num", 1000);
			iter_count = Json::Value(j, "first_iter", 0);
			output_base_dir = Json::Value(j, "output_base_dir", std::string("output"));
			base_path = bf::path(output_base_dir);
		}
	};
}