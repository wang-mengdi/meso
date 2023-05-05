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
		if (base_path.is_relative()) {
			base_path = bf::current_path() / base_path;
		}

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

	bf::path DriverMetaData::Current_VTS_Path(const std::string identifier)
	{
		return base_path / bf::path(fmt::format("{}{:04d}.vts", identifier, current_frame));
	}

	bf::path DriverMetaData::Current_OBJ_Path(const std::string identifier)
	{
		return base_path / bf::path(fmt::format("{}{:04d}.obj", identifier, current_frame));
	}

	bool DriverMetaData::Should_Snapshot(void)
	{
		return current_frame != 0 && snapshot_stride != 0 && current_frame % snapshot_stride == 0;
	}

	bf::path DriverMetaData::Snapshot_Base_Path(void)
	{
		return base_path / bf::path("snapshots");
	}

	bf::path DriverMetaData::Snapshot_Path(int frame)
	{
		return Snapshot_Base_Path() / bf::path(fmt::format("{:04d}", frame));
	}

	bf::path DriverMetaData::Current_Snapshot_Path(void)
	{
		return Snapshot_Path(current_frame);
	}

	int DriverMetaData::Last_Snapshot_Frame(int start_frame)
	{
		bf::path snap_base = Snapshot_Base_Path();
		if (!bf::is_directory(snap_base)) return 0;//no snapshots are there

		std::vector<int> snapshots;
		for (bf::directory_iterator itr(snap_base); itr != bf::directory_iterator(); ++itr) {
			if (bf::is_directory(itr->status())) {
				std::string filename = itr->path().filename().stem().string();
				snapshots.push_back(std::stoi(filename));
			}
		}

		std::sort(snapshots.begin(), snapshots.end());
		//find the first element >= start_frame
		auto it = std::lower_bound(snapshots.begin(), snapshots.end(), start_frame);
		if (it != snapshots.begin()) {
			it--;
			return *it;
		}
		else {
			return 0;//no snapshot is read
		}
	}
}