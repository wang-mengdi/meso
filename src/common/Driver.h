//////////////////////////////////////////////////////////////////////////
// Simulator Driver
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Common.h"
#include "AuxFunc.h"
#include "Json.h"

//#include <iostream>
//#include <filesystem>

//namespace fs = std::filesystem;

namespace Meso {

	template<class Simulator>
	class Driver {
	public:
		Simulator simulator;

		int fps = 25;
		real cfl = 1.0;
		real time_per_frame = 0.04;
		real min_step_frame_fraction = 0;//if set to 0.1, it means the minimal iteration time is 0.1*time_per_frame

		int first_frame;
		int last_frame;

		std::string output_base_dir;
		
		void Init(json& j) {
			fps = Json::Value(j, "fps", 25);
			cfl = Json::Value(j, "cfl", 1.0);
			time_per_frame = 1.0 / fps;
			min_step_frame_fraction = Json::Value(j, "min_step_frame_fraction", (real)0);
			first_frame = Json::Value(j, "first_frame", 0);
			last_frame = Json::Value(j, "last_frame", fps * 10);
			output_base_dir = Json::Value(j, "output_base_dir", "output");
		}
		void Time_At_Frame(const int frame) {
			return frame * time_per_frame;
		}
		//at the beginning the system is at the status of start_frame
		//will output all frames in [start_frame, end_frame]
		void Advance(int start_frame, int end_frame) {
			simulator.Output(output_base_dir, start_frame);
			for (int current_frame = start_frame; current_frame < end_frame; current_frame++) {
				int next_frame = current_frame + 1;
				real current_time = Time_At_Frame(current_frame);
				real next_time = Time_At_Frame(next_frame);
				while (true) {
					//can return an inf
					real dt = simulator.CFL_Time(cfl);
					dt = MathFunc::Clamp(dt, min_step_time, time_per_frame);
					if (current_time + dt >= next_time) {
						dt = next_time - current_time;
						simulator.Advance(current_frame, current_time, dt);
						break;
					}
					else simulator.Advance(current_frame, current_time, dt);
				}
				simulator.Output(output_base_dir, next_frame);
			}
		}

		template<class Initializer>
		void Run(json& j, Initializer& initializer) {
			Init(j.at("driver"));
			initializer.Init(j.at("initializer"));
			initializer.Apply(simulator);
			Advance(first_frame, last_frame);
		}
	};
}