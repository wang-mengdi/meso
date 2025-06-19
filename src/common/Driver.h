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
		void Print_Frame_Info(Timer& frame_timer, const DriverMetaData& meta_data);
		//will change timer
		void Print_Iteration_Info(Timer& iteration_timer, const real dt, const real running_cfl, const real current_time, const real frame_time);

		//simulate from first_frame to last_frame
		//output frames [first_frame, last_frame]
		void Advance(Simulator& simulator, DriverMetaData& meta_data);

		void Run(json& j, Simulator& simulator);

		template<class Initializer, class TSimulator>
		void Initialize_And_Run(json& j, Initializer& scene, TSimulator& simulator) {
			Info("Driver::Initialize_And_Run parse json: \n{}", j.dump(2));
			DriverMetaData meta_data;
			meta_data.Init(j.at("driver"));
			scene.Apply(j, simulator);
			FileFunc::Create_Directory(meta_data.output_base_dir);
			fs::path dump_file = fs::path(meta_data.output_base_dir) / fs::path("config.json");
			std::ofstream dump_output(dump_file.string());
			dump_output <<std::setw(4)<< j;
			dump_output.close();
			Advance(simulator, meta_data);
		}
	};
}