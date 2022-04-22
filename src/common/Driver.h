//////////////////////////////////////////////////////////////////////////
// Simulator Driver
// Copyright (c) (2022-), Bo Zhu, Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once

#include "Common.h"
#include "AuxFunc.h"
#include "Json.h"
#include "Simulator.h"

#include <boost/filesystem.hpp>
namespace bf = boost::filesystem;

namespace Meso {

	class Driver {
	public:
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
		real Time_At_Frame(const int frame) {
			return frame * time_per_frame;
		}
		//at the beginning the system is at the status of start_frame
		//will output all frames in [start_frame, end_frame]
		void Advance(Simulator &simulator, int start_frame, int end_frame) {
			bf::path base_path(output_base_dir);
			bf::path frame_dir(std::to_string(start_frame));
			simulator.Output(base_path.string(), (base_path / frame_dir).string());
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
				frame_dir = bf::path(std::to_string(next_frame));
				simulator.Output(base_path.string(), (base_path / frame_dir).string());
			}
		}

		template<class Initializer>
		void Run(json& j, Initializer& scene, Simulator& simulator) {
			Init(j.at("driver"));
			scene.Apply(j.at("scene"), simulator);
			Advance(simulator, first_frame, last_frame);
		}
	};
}