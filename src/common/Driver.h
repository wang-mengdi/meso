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
#include "MetaData.h"
#include <fstream>

namespace Meso {

	class Driver {
	public:
		//will change timer
		void Print_Frame_Info(Timer &frame_timer, const DriverMetaData& meta_data) {
			int total_frames = meta_data.last_frame - meta_data.first_frame;
			int done_frames = meta_data.current_frame - meta_data.first_frame;
			real frame_seconds = frame_timer.Lap_Time();
			real completed_seconds = frame_timer.Total_Time();
			real eta = completed_seconds * (total_frames - done_frames) / done_frames;
			Info("Frame {} in {}-{} done in {:.3f}s, ETA {:.3f}/{:.3f}s", meta_data.current_frame, meta_data.first_frame, meta_data.last_frame, frame_seconds, eta, completed_seconds + eta);
		}
		//will change timer
		void Print_Iteration_Info(Timer& iteration_timer, const real dt, const real current_time, const real frame_time) {
			real step_seconds = iteration_timer.Lap_Time();
			real completed_seconds = iteration_timer.Total_Time();
			Info("Iteration {:.5f}/{:.5f}s, cost {:.3f}s, remaining {:.3f}s", dt, frame_time, step_seconds, completed_seconds * (frame_time - current_time) / current_time);
		}

		//at the beginning the system is at the status of start_frame
		//will output all frames in [start_frame, end_frame]
		void Advance(Simulator &simulator, DriverMetaData& meta_data) {
			Timer frame_timer;
			bf::path base_path(meta_data.output_base_dir);
			FileFunc::Create_Directory(base_path);

			//output first frame
			Print_Frame_Info(frame_timer, meta_data);
			Info("Output frame {} to {}", meta_data.current_frame, meta_data.output_base_dir);
			simulator.Output(meta_data);
			meta_data.current_frame += 1;

			//run frames
			for (int& current_frame = meta_data.current_frame; current_frame <= meta_data.last_frame; current_frame++) {
				Timer iter_timer;
				int next_frame = current_frame + 1;
				meta_data.current_time = meta_data.Time_At_Frame(current_frame);
				real frame_start_time = meta_data.current_time;
				real next_time = meta_data.Time_At_Frame(next_frame);
				while (true) {
					//can return an inf
					meta_data.dt = simulator.CFL_Time(meta_data.cfl);
					meta_data.dt = MathFunc::Clamp(meta_data.dt, meta_data.min_step_frame_fraction * meta_data.time_per_frame, meta_data.time_per_frame);
					bool last_iter = false;
					if (meta_data.current_time + meta_data.dt >= next_time) {
						meta_data.dt = next_time - meta_data.current_time;
						last_iter = true;
					}
					else if (meta_data.current_time + 2 * meta_data.dt >= next_time) {
						meta_data.dt = (real).5 * (next_time - meta_data.current_time);
					}


					simulator.Advance(meta_data);
					meta_data.current_time += meta_data.dt;
					Print_Iteration_Info(iter_timer, meta_data.dt, meta_data.current_time - frame_start_time, meta_data.time_per_frame);
					if (last_iter) break;
				}

				//output current frame
				Print_Frame_Info(frame_timer, meta_data);
				Info("Output frame {} to {}", meta_data.current_frame, meta_data.output_base_dir);
				simulator.Output(meta_data);
			}
		}

		template<class Initializer, class TSimulator>
		void Initialize_And_Run(json& j, Initializer& scene, TSimulator& simulator) {
			Info("Driver::Run parse json: \n{}", j.dump(2));
			DriverMetaData meta_data;
			meta_data.Init(j.at("driver"));
			scene.Apply(j.at("scene"), simulator);
			FileFunc::Create_Directory(meta_data.output_base_dir);
			bf::path dump_file = bf::path(meta_data.output_base_dir) / bf::path("config.json");
			std::ofstream dump_output(dump_file.string());
			dump_output <<std::setw(4)<< j;
			dump_output.close();
			Advance(simulator, meta_data);
		}
	};
}