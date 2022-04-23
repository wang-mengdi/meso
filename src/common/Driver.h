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
#include "Timer.h"

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
			output_base_dir = Json::Value(j, "output_base_dir", std::string("output"));
		}
		real Time_At_Frame(const int frame) {
			return frame * time_per_frame;
		}
		//will change timer
		void Print_Frame_Info(Timer &frame_timer, const int frame, const int start_frame, const int end_frame) {
			int total_frames = end_frame - start_frame;
			int done_frames = frame - start_frame;
			real frame_seconds = frame_timer.Lap_Time();
			real completed_seconds = frame_timer.Total_Time();
			real eta = completed_seconds * total_frames / done_frames;
			Info("Frame {} in {}-{} done in {:.3f}s, ETA {:.3f}/{:.3f}s", frame, start_frame, end_frame, frame_seconds, eta, completed_seconds + eta);
		}
		//will change timer
		void Print_Iteration_Info(Timer& iteration_timer, const real dt, const real current_time, const real frame_time) {
			real step_seconds = iteration_timer.Lap_Time();
			real completed_seconds = iteration_timer.Total_Time();
			Info("Iteration {:.5f}/{:.5f}s, cost {:.3f}s, ETA {:.3f}s", dt, frame_time, step_seconds, completed_seconds * frame_time / current_time);
		}

		void Output(Simulator &simulator, bf::path base_path, const int frame) {
			//bf::path frame_dir(std::to_string(frame));
			//bf::path frame_path = base_path / frame_dir;
			//IOFunc::Create_Directory(frame_path);
			Info("Output frame {} to {}", frame, base_path.string());
			simulator.Output(base_path, frame);
		}

		//at the beginning the system is at the status of start_frame
		//will output all frames in [start_frame, end_frame]
		void Advance(Simulator &simulator, int start_frame, int end_frame) {
			Timer frame_timer;
			bf::path base_path(output_base_dir);
			IOFunc::Create_Directory(base_path);

			Print_Frame_Info(frame_timer, start_frame, start_frame, end_frame);
			Output(simulator, base_path, start_frame);
			for (int current_frame = start_frame; current_frame <= end_frame; current_frame++) {
				Timer iter_timer;
				int next_frame = current_frame + 1;
				real current_time = Time_At_Frame(current_frame);
				real frame_start_time = current_time;
				real next_time = Time_At_Frame(next_frame);
				while (true) {
					//can return an inf
					real dt = simulator.CFL_Time(cfl);
					dt = MathFunc::Clamp(dt, min_step_frame_fraction * time_per_frame, time_per_frame);
					bool last_iter = false;
					if (current_time + dt >= next_time) {
						dt = next_time - current_time;
						last_iter = true;
					}
					simulator.Advance(current_frame, current_time, dt);
					current_time += dt;
					Print_Iteration_Info(iter_timer, dt, current_time - frame_start_time, time_per_frame);
					if (last_iter) break;
				}
				Print_Frame_Info(frame_timer, current_frame, start_frame, end_frame);
				Output(simulator, base_path, current_frame);
			}
		}

		template<class Initializer, class TSimulator>
		void Run(json& j, Initializer& scene, TSimulator& simulator) {
			Info("Driver::Run parse json: \n{}", j.dump(2));
			Init(j.at("driver"));
			scene.Apply(j.at("scene"), simulator);
			Advance(simulator, first_frame, last_frame);
		}
	};
}