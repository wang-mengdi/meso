//////////////////////////////////////////////////////////////////////////
// Metadata
// Copyright (c) (2022-), Mengdi Wang, Yuchen Sun
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
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

		~DriverMetaData();
		real Time_At_Frame(int frame);
		void Init(json& j);
		void Append_Output_Thread(std::shared_ptr<std::thread> thread_ptr);

		//path of an output .vts file at current frame, for example pressure, velocity
		bf::path Current_VTS_Path(const std::string identifier);
		bf::path Current_VTU_Path(const std::string identifier);
		bf::path Current_OBJ_Path(const std::string identifier);

		//snapshot things
		bool Should_Snapshot(void);
		bf::path Snapshot_Base_Path(void);
		bf::path Snapshot_Path(int frame);
		bf::path Current_Snapshot_Path(void);//snapshot path of current frame
		int Last_Snapshot_Frame(int start_frame);//last snapshotted frame < start_frame
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