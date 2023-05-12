//////////////////////////////////////////////////////////////////////////
// Timer
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include <chrono>
#include <ctime>
#include "Common.h"
#include "Constants.h"
#include <cuda_runtime.h>

namespace Meso {
	class TimerRecord {
	public:
		double total;//total time
		double cur;//time of the current loop
		int num;//number of loops
		std::string name;
		TimerRecord(std::string _name) :total(0), cur(0), num(0), name(_name) {}
		void Add(double cur_time);//add a loop
		std::pair<double, double> Profile(void);//first:current time. second: avg time
	};

	class Timer {
	public:
		std::chrono::time_point<std::chrono::system_clock> total_start;
		std::chrono::time_point<std::chrono::system_clock> lap_start;
		std::chrono::time_point<std::chrono::system_clock> loop_start;
		std::map<std::string, int> name_index;
		Array<TimerRecord> records;

		Timer() {
			Reset();
		}
		void Reset(void);
		//total time
		real Total_Time(const real unit = PhysicalUnits::s);
		//lap time, and reset the lap clock
		real Lap_Time(const real unit = PhysicalUnits::s);
		//loop functions
		void Begin_Loop(void);
		void Record(const std::string& name, const real unit = PhysicalUnits::s);
		void Output_Profile(std::ostream& out = std::cout);
	};

	class GPUTimer {
	public:
		cudaEvent_t start, stop;

		GPUTimer() {
			cudaEventCreate(&start);
			cudaEventCreate(&stop);
		}
		~GPUTimer() {
			cudaEventDestroy(start);
			cudaEventDestroy(stop);
		}

		void Start(void);
		void Stop(const std::string& output_message);
	};
}