#pragma once
#include "Common.h"
#include "Json.h"
#include <fstream>
#include "Timer.h"
#include <queue>
#include <thread>
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
		bool output_each_step;				//whether or not output each step, ignoring the frame rate, each frame has duration cfl_time
		real total_time;					//total time duration, only used when output_each_step
		int first_frame;
		int last_frame;
		int snapshot_stride = 0;
		
		int output_queue_size = 10;

		//queue of threads for output
		std::queue<std::shared_ptr<std::thread>> output_threads;

		//fill for every time step
		int current_frame;
		real current_time;
		real dt;
		real running_cfl;//actual cfl number of simulation

		~DriverMetaData() {
			while (!output_threads.empty()) {
				auto join_ptr = output_threads.front();
				output_threads.pop();
				join_ptr->join();
			}
		}

		real Time_At_Frame(int frame) {
			return frame * time_per_frame;
		}

		void Init(json& j) {
			output_base_dir = Json::Value(j, "output_base_dir", std::string("output"));
			base_path = bf::path(output_base_dir);

			output_each_step = Json::Value(j, "output_each_step", false);
			if (output_each_step) {
				total_time = Json::Value(j, "total_time", (real)1);
			}
			else {
				fps = Json::Value(j, "fps", 25);
				cfl = Json::Value(j, "cfl", (real)1.0);
				time_per_frame = 1.0 / fps;
				min_step_frame_fraction = Json::Value(j, "min_step_frame_fraction", (real)0);
				last_frame = Json::Value(j, "last_frame", fps * 10);
				dt = 1.0 / fps;
			}

			first_frame = Json::Value(j, "first_frame", 0);
			snapshot_stride = Json::Value(j, "snapshot_stride", 0);
			
			output_queue_size = Json::Value(j, "queue_size", 10);

			current_frame = first_frame;
			current_time = Time_At_Frame(current_frame);
		}

		void Append_Output_Thread(std::shared_ptr<std::thread> thread_ptr) {
			while (output_threads.size() > output_queue_size) {
				auto join_ptr = output_threads.front(); 
				output_threads.pop();
				join_ptr->join();
			}
			output_threads.push(thread_ptr);
		}
	};

	class OptimizerDriverMetaData : public MetaData
	{
	public:
		real tol;
		int iter_count;
		int max_iter_num;
		bool verbose;
		Timer timer;
		void Init(json& j) {
			tol = Json::Value(j, "tol", (real) 1e-7);
			max_iter_num = Json::Value(j, "max_iter_num", 1000);
			iter_count = Json::Value(j, "first_iter", 0);
			output_base_dir = Json::Value(j, "output_base_dir", std::string("output"));
			verbose = Json::Value(j, "verbose", true);
			base_path = bf::path(output_base_dir);
		}
	};
}