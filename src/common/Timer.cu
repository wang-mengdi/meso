#include "Timer.h"

namespace Meso {
	//time record
	void TimerRecord::Add(double cur_time)
	{
		total += cur_time;
		cur = cur_time;
		num++;
	}

	std::pair<double, double> TimerRecord::Profile(void)
	{
		return std::pair<double, double>(cur, total / num);
	}

	//timer
	void Timer::Reset(void)
	{
		total_start = std::chrono::system_clock::now();
		lap_start = std::chrono::system_clock::now();
	}

	real Timer::Total_Time(const real unit)
	{
		std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
		std::chrono::duration<double, std::ratio<1> > elapse = end - total_start;
		return elapse.count() / unit;
	}

	real Timer::Lap_Time(const real unit)
	{
		std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
		std::chrono::duration<double, std::ratio<1> > elapse = end - lap_start;
		lap_start = std::chrono::system_clock::now();

		return elapse.count() / unit;
	}

	void Timer::Begin_Loop(void)
	{
		loop_start = std::chrono::system_clock::now();
		Reset();
	}

	void Timer::Record(const std::string& name, const real unit)
	{
		int k = -1;
		auto iter = name_index.find(name);
		if (iter == name_index.end()) {
			k = records.size();
			name_index[name] = k;
			records.push_back(TimerRecord(name));//push a zero record
		}
		else {
			k = iter->second;
		}
		records[k].Add(Lap_Time(unit));
	}

	void Timer::Output_Profile(std::ostream& out)
	{
		std::string str = "#     Time record cur(avg): ";
		for (int i = 0; i < records.size(); i++) {
			std::pair<double, double> rec = records[i].Profile();
			str += fmt::format("{}: {:.3f}({:.3f}), ", records[i].name, rec.first, rec.second);
		}
		auto loop_end = std::chrono::system_clock::now();
		std::chrono::duration<double, std::ratio<1> > loop_elapse = loop_end - loop_start;//in seconds
		double loop_time_sec = loop_elapse.count() / PhysicalUnits::s;
		str += fmt::format("total: {:.3f}", loop_time_sec);
		out << str << "\n";
	}

}