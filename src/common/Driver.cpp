#include "Driver.h"

namespace Meso {
	void Driver::Print_Frame_Info(Timer& frame_timer, const DriverMetaData& meta_data) {
		int total_frames = meta_data.last_frame - meta_data.first_frame;
		int done_frames = meta_data.current_frame - meta_data.first_frame;
		real frame_seconds = frame_timer.Lap_Time();
		real completed_seconds = frame_timer.Total_Time();
		real eta = completed_seconds * (total_frames - done_frames) / done_frames;
		Info("Frame {} in {}-{} done in {:.3f}s, ETA {:.3f}/{:.3f}s", meta_data.current_frame, meta_data.first_frame, meta_data.last_frame, frame_seconds, eta, completed_seconds + eta);
	}

	void Driver::Print_Iteration_Info(Timer& iteration_timer, const real dt, const real running_cfl, const real current_time, const real frame_time)
	{
		real step_seconds = iteration_timer.Lap_Time();
		real completed_seconds = iteration_timer.Total_Time();
		Info("Iteration {:.5f}/{:.5f}s with CFL={:.3f}, cost {:.3f}s, remaining {:.3f}s", dt, frame_time, running_cfl, step_seconds, completed_seconds * (frame_time - current_time) / current_time);
	}

	void Driver::Advance(Simulator& simulator, DriverMetaData& meta_data) {
		Timer frame_timer;
		FileFunc::Create_Directory(meta_data.base_path);

		//try to load snapshot
		if (meta_data.first_frame != 0) {
			int last_snapshot = meta_data.Last_Snapshot_Frame(meta_data.first_frame);
			if (last_snapshot != 0) {
				Info("Found snapshot at frame {}, load snapshot from {} and run from frame {}", last_snapshot, meta_data.Snapshot_Path(last_snapshot).string(), last_snapshot + 1);
				meta_data.current_frame = last_snapshot;
				//load the frame before first frame
				simulator.Load_Frame(meta_data);
				meta_data.current_frame++;
			}
			else {
				Info("No snapshots are found, automaitcally run from frame 0");
				meta_data.first_frame = 0;
			}
		}

		//run from 0
		if (meta_data.first_frame == 0) {
			//output first frame
			meta_data.current_frame = 0;
			Print_Frame_Info(frame_timer, meta_data);
			Info("Output frame {} to {}", meta_data.current_frame, meta_data.base_path.string());
			simulator.Output(meta_data);
			meta_data.current_frame = 1;
		}

		//run frames
		for (int& current_frame = meta_data.current_frame; current_frame <= meta_data.last_frame; current_frame++) {
			Timer iter_timer;
			int next_frame = current_frame + 1;
			meta_data.current_time = meta_data.Time_At_Frame(current_frame);
			real frame_start_time = meta_data.current_time;
			real next_time = meta_data.Time_At_Frame(next_frame);
			while (true) {
				//can return an inf
				real math_dt = simulator.CFL_Time(meta_data.cfl);
				meta_data.dt = MathFunc::Clamp(math_dt, meta_data.min_step_frame_fraction * meta_data.time_per_frame, meta_data.time_per_frame);
				meta_data.running_cfl = meta_data.cfl * meta_data.dt / math_dt;
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
				Print_Iteration_Info(iter_timer, meta_data.dt, meta_data.running_cfl, meta_data.current_time - frame_start_time, meta_data.time_per_frame);
				if (last_iter) break;
			}

			//output current frame
			Print_Frame_Info(frame_timer, meta_data);
			Info("Output frame {} to {}", meta_data.current_frame, meta_data.base_path.string());
			simulator.Output(meta_data);
			Info("==========================================");
		}
	}

	void Driver::Run(json& j, Simulator& simulator) {
		Info("Driver::Run parse json: \n{}", j.dump(2));
		DriverMetaData meta_data;
		meta_data.Init(j.at("driver"));
		FileFunc::Create_Directory(meta_data.base_path.string());
		fs::path dump_file = meta_data.base_path / fs::path("config.json");
		std::ofstream dump_output(dump_file.string());
		dump_output << std::setw(4) << j;
		dump_output.close();
		Advance(simulator, meta_data);
	}
}