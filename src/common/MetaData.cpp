#include "MetaData.h"

namespace Meso {
	DriverMetaData::~DriverMetaData() {
		while (!output_threads.empty()) {
			auto join_ptr = output_threads.front();
			output_threads.pop();
			join_ptr->join();
		}
	}

	real DriverMetaData::Time_At_Frame(int frame) {
		return frame * time_per_frame;
	}

	void DriverMetaData::Init(json& j) {
		output_base_dir = Json::Value(j, "output_base_dir", std::string("output"));
		base_path = bf::path(output_base_dir);

		output_each_step = Json::Value(j, "output_each_step", false);
		if (output_each_step) {
			total_time = Json::Value(j, "total_time", (real)1);
		}
		else {
			fps = Json::Value(j, "fps", 25);
			cfl = Json::Value(j, "cfl", (real)1.0);
			time_per_frame = 1.0 / fps;
			min_step_frame_fraction = Json::Value(j, "min_step_frame_fraction", (real)0);
			last_frame = Json::Value(j, "last_frame", fps * 10);
			dt = 1.0 / fps;
		}

		first_frame = Json::Value(j, "first_frame", 0);
		snapshot_stride = Json::Value(j, "snapshot_stride", 0);

		output_queue_size = Json::Value(j, "queue_size", 10);

		current_frame = first_frame;
		current_time = Time_At_Frame(current_frame);
	}

	void DriverMetaData::Append_Output_Thread(std::shared_ptr<std::thread> thread_ptr) {
		while (output_threads.size() > output_queue_size) {
			auto join_ptr = output_threads.front();
			output_threads.pop();
			join_ptr->join();
		}
		output_threads.push(thread_ptr);
	}
}