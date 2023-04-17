#include "IOFunc.h"
namespace Meso {
	namespace VTKFunc {
		void Write_VTS(vtkStructuredGrid& vtk_grid, std::string file_name) {
			// setup VTK
			vtkNew<vtkXMLStructuredGridWriter> writer;

#if (VTK_MAJOR_VERSION >=6)
			writer->SetInputData(&vtk_grid);
#else
			writer->SetInput(&vtk_grid);
#endif

			writer->SetFileName(file_name.c_str());
			writer->SetDataModeToBinary();
			writer->Write();
		}

		void Write_VTP(vtkPolyData& vtp, std::string file_name) {
			if (vtp.GetNumberOfPoints() == 0) { Pass("here no points"); return; }
			// setup VTK
			vtkNew<vtkXMLPolyDataWriter> writer;

			writer->SetFileName(file_name.c_str());
#if (VTK_MAJOR_VERSION >=6)
			writer->SetInputData(&vtp);
#else
			writer->SetInput(&vtk_grid);
#endif
			writer->SetDataModeToBinary();
			writer->Write();
		}
	}
}