#include "Common.h"

namespace Meso {

	void Info(const std::string& str)
	{
		Info(str.c_str());
        //Info(str);
	}

	void Warn(const std::string& str)
	{
		Warn(str.c_str());
	}

    void Check_Cuda_Memory(const std::string str) {
        size_t free_byte;

        size_t total_byte;

        auto cuda_status = cudaMemGetInfo(&free_byte, &total_byte);

        if (cudaSuccess != cuda_status) {

            Error("{} cudaMemGetInfo fails: {} ", str, cudaGetErrorString(cuda_status));

            exit(1);

        }

        double free_db = (double)free_byte;

        double total_db = (double)total_byte;

        double used_db = total_db - free_db;

        Info("{} GPU memory usage: used = {}, free = {} MB, total = {} MB\n",

            str, used_db / 1024.0 / 1024.0, free_db / 1024.0 / 1024.0, total_db / 1024.0 / 1024.0);
    }

}