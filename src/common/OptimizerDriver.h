#pragma once 
#include "common.h"
#include "Json.h"
#include "Timer.h"
#include "MetaData.h"
#include "Optimizer.h"
#include "AuxFunc.h"

namespace Meso{
	class OptimizerDriver
	{
	public:
		OptimizerDriver(){}

		//timer changed
		void Print_Iteration_Info(Timer& iteration_timer,const OptimizerDriverMetaData& meta_data) {
			real step_seconds = iteration_timer.Lap_Time();
			real completed_seconds = iteration_timer.Total_Time();
			Info("Iteration {} cost {:.3f}s, {}s passed", meta_data.iter_count, step_seconds, completed_seconds);
		}

		void Output(Optimizer& optimizer, OptimizerDriverMetaData& meta_data) {
			Info("Output iter {} to {}", meta_data.iter_count, meta_data.base_path.string());
			optimizer.Output(meta_data);
			meta_data.data_output.flush();
		}

		void Advance(Optimizer& optimizer, OptimizerDriverMetaData& meta_data) {
			Timer iter_timer;
			fs::path base_path(meta_data.output_base_dir);
			FileFunc::Create_Directory(base_path);

			Output(optimizer, meta_data);
			meta_data.iter_count++;
			meta_data.timer.Begin_Loop();
			while (!optimizer.Is_Converged(meta_data) && meta_data.iter_count <= meta_data.max_iter_num) {
				optimizer.Optimize(meta_data);
				if (meta_data.verbose) { Print_Iteration_Info(iter_timer, meta_data); }
				Output(optimizer, meta_data);
				if (meta_data.verbose) { meta_data.timer.Output_Profile(); }
				meta_data.iter_count++;
			}
			fs::path timing_file = fs::path(meta_data.output_base_dir) / fs::path("timing.txt");
			std::ofstream timing_output(timing_file.string());
			meta_data.timer.Output_Profile(timing_output);
			timing_output.close();
			meta_data.data_output.close();
			if (meta_data.iter_count== meta_data.max_iter_num) {
				Error("Optimizer doesn't converge after {} iters!", meta_data.iter_count);
			}
			else {
				Info("Optimizer converges after {} iters", meta_data.iter_count - 1);
			}
		}

		template<class Initializer, class TOptimizer>
		void Initialize_And_Run(json& j, Initializer& scene, TOptimizer& optimizer) {
			OptimizerDriverMetaData meta_data=Initialize(j, scene, optimizer);
			Advance(optimizer, meta_data);
		}

		template<class Initializer, class TOptimizer>
		OptimizerDriverMetaData Initialize(json& j, Initializer& scene, TOptimizer& optimizer) {
			Info("OptimizerDriver::Run parse json: \n{}", j.dump(2));
			OptimizerDriverMetaData meta_data;
			meta_data.Init(j.at("optimizer"));
			scene.Apply(j.at("scene"), optimizer);
			FileFunc::Create_Directory(meta_data.output_base_dir);
			fs::path dump_file = fs::path(meta_data.output_base_dir) / fs::path("config.json");
			std::ofstream dump_output(dump_file.string());
			dump_output << std::setw(4) << j;
			dump_output.close();
			fs::path output_data_path = fs::path(meta_data.output_base_dir) / fs::path("output_data.csv");
			meta_data.data_output.open(output_data_path.string());
			return meta_data;
		}
	};
}