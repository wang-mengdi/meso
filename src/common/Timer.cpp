#include "Timer.h"

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
