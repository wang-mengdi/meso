#pragma once
#include "Common.h"
#include "Json.h"

namespace Meso {
	class MetaData
	{
	public:
		std::string output_base_dir;
	};

	class DriverMetaData : public MetaData
	{
	public:
		int fps = 25;
		real cfl = 1.0;
		real time_per_frame = 0.04;
		real min_step_frame_fraction = 0;	//if set to 0.1, it means the minimal iteration time is 0.1*time_per_frame

		int frame;							//current frame
		int first_frame;
		int last_frame;

		real Time_At_Frame(int frame) {
			return frame * time_per_frame;
		}

		void Init(json& j) {
			fps = Json::Value(j, "fps", 25);
			cfl = Json::Value(j, "cfl", 1.0);
			time_per_frame = 1.0 / fps;
			min_step_frame_fraction = Json::Value(j, "min_step_frame_fraction", (real)0);
			first_frame = Json::Value(j, "first_frame", 0);
			frame = first_frame;
			last_frame = Json::Value(j, "last_frame", fps * 10);
			output_base_dir = Json::Value(j, "output_base_dir", std::string("output"));
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
		}
	};
}