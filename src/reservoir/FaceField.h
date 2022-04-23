//////////////////////////////////////////////////////////////////////////
// Face data on MacGrid
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Grid.h"
#include "Interpolation.h"
#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>


namespace Meso {

	template<class T, int d, DataHolder side = DataHolder::HOST>
	class FaceField {
		Typedef_VectorD(d);
	public:
		Grid<d> grid;
		//std::array<Array<T, side>, d> face_data;
		std::shared_ptr<Array<T, side>> face_data[3] = { nullptr,nullptr,nullptr };
		FaceField() {}
		FaceField(const Grid<d>& _grid) { Init(_grid); }
		FaceField(const Grid<d>& _grid, const T value) { Init(_grid);  Fill(value); }
		template<DataHolder side1> FaceField(const FaceField<T, d, side1>& f1) { Init(f1); }
		void Fill(const T value) { for (int axis = 0; axis < d; axis++) ArrayFunc::Fill(*face_data[axis], value); }
		void Init(const Grid<d>& _grid, const T value) { Init(_grid); Fill(value); }
		void Init(const Grid<d>& _grid) {
			grid = _grid;
			for (int axis = 0; axis < d; axis++) {
				int n = grid.Face_DoF(axis);
				if (face_data[axis] == nullptr) face_data[axis] = std::make_shared<Array<T, side>>(n);
				else face_data[axis]->resize(n);
				checkCudaErrors(cudaGetLastError());
			}
		}
		template<DataHolder side1>
		void Init(const FaceField<T, d, side1>& f1) {
			Deep_Copy(f1);
		}

		template<DataHolder side1> 
		void Deep_Copy(const FaceField<T, d, side1>& f1) {
			Init(f1.grid);
			for (int i = 0; i < d; i++) {
				//deep copy
				*face_data[i] = f1.Data(i);
			}
		}

		template<DataHolder side1>
		FaceField<T, d, side>& operator = (const FaceField<T, d, side1>& f1) {
			Deep_Copy(f1);
			return *this;
		}
		template<DataHolder side1>
		FaceField<T, d, side>& operator = (FaceField<T, d, side1>& f1) {
			Deep_Copy(f1);
			return *this;
		}
		inline T& operator()(const int axis, const VectorDi face) { return (*(face_data[axis]))[grid.Face_Index(axis, face)]; }
		inline const T& operator()(int axis, const VectorDi face) const { return (*(face_data[axis]))[grid.Face_Index(axis, face)]; }
		inline const T Get(int axis, const VectorDi face) const { return (*(face_data[axis]))[grid.Face_Index(axis, face)]; }
		void operator += (const FaceField<T, d, side>& f1) {
			for (int axis = 0; axis < d; axis++) {
				ArrayFunc::Add(Data(axis), f1.Data(axis));
			}
		}

		constexpr Array<T, side>& Data(const int axis)noexcept { return *face_data[axis]; }
		constexpr const Array<T, side>& Data(const int axis)const noexcept { return *face_data[axis]; }
		constexpr T* Data_Ptr(const int axis) noexcept { return face_data[axis] == nullptr ? nullptr : thrust::raw_pointer_cast(face_data[axis]->data()); }
		constexpr const T* Data_Ptr(const int axis) const noexcept {
			return face_data[axis] == nullptr ? nullptr : thrust::raw_pointer_cast(face_data[axis]->data());
		}

		void operator *= (const FaceField<T, d, side>& f1) {
			for (int axis = 0; axis < d; axis++) {
				ArrayFunc::Multiply(Data(axis), f1.Data(axis));
			}
		}
		template<class T1> void operator *= (const T1 a) {
			for (int axis = 0; axis < d; axis++) ArrayFunc::Multiply_Scalar(Data(axis), a);
		}

		T Max_Abs(void) {
			real max_val = 0;
			for (int axis = 0; axis < d; axis++) {
				max_val = std::max<T>(max_val, ArrayFunc::Max_Abs<T>(Data(axis)));
			}
			return max_val;
		}

		void Output_Vtk(const std::string file_name) {
			// Reading data from the file. ">>" is a operator of "ifstream"
			int nx, ny, nz;
			if constexpr (d == 2) {
				nx = grid.counts[0];
				ny = grid.counts[1];
				nz = 1;
			}
			else {
				nx = grid.counts[0];
				ny = grid.counts[1];
				nz = grid.counts[2];
			}

			// setup VTK
			vtkXMLStructuredGridWriter* writer = vtkXMLStructuredGridWriter::New();
			vtkStructuredGrid* structured_grid = vtkStructuredGrid::New();
			structured_grid->SetDimensions(nx, ny, nz);
			vtkPoints* nodes = vtkPoints::New();
			nodes->Allocate(nx * ny * nz);
			//vtkDoubleArray* prsArray = vtkDoubleArray::New();
			//prsArray->SetNumberOfComponents(1);
			//prsArray->SetName("Pressure");
			vtkDoubleArray* velArray = vtkDoubleArray::New();
			velArray->SetName("Velocity");
			velArray->SetNumberOfComponents(3);

			FaceField<T, d> field_host = *this;

			for (int k = 0; k < nz; k++) {
				for (int j = 0; j < ny; j++) {
					for (int i = 0; i < nx; i++) {
						real x, y, z, u, v, w, p;
						//ifile >> x >> y >> z >> u >> v >> w >> p;
						VectorDi cell = VectorFunc::Vi<d>(i, j, k);
						//VectorD pos = VectorFunc::V<d>(i, j, k);
						VectorD pos = grid.Position(cell);
						Vector3 pos3 = VectorFunc::V<3>(pos);

						//nodes->InsertNextPoint(x, y, z);
						nodes->InsertNextPoint(pos3[0], pos3[1], pos3[2]);
						VectorD vel = IntpLinear::Face_Vector<T, d, HOST>(field_host, pos);
						Vector3 vel3 = VectorFunc::V<3>(vel);
						//velArray->InsertNextTuple3(u, v, w);
						//velArray->InsertNextTuple3(vel3[0], vel3[1], vel3[2]);
						velArray->InsertNextTuple3(vel3[0], vel3[1], vel3[2]);
					}
				}
			}

			structured_grid->SetPoints(nodes);
			structured_grid->GetPointData()->AddArray(velArray);
			structured_grid->GetPointData()->SetActiveVectors("velocity");
			//structured_grid->GetPointData()->AddArray(prsArray);

#if (VTK_MAJOR_VERSION >=6)
			writer->SetInputData(structured_grid);
#else
			writer->SetInput(structured_grid);
#endif

			writer->SetFileName(file_name.c_str());
			writer->SetDataModeToAscii();
			writer->Write();

			structured_grid->Delete();
			writer->Delete();
			nodes->Delete();
			velArray->Delete();
			//prsArray->Delete();
		}

		template<class IFFunc>
		void Iterate_Faces(IFFunc f) {
			for (int axis = 0; axis < d; axis++) {
				int n = grid.Face_DoF(axis);
				for (int i = 0; i < n; i++) {
					VectorDi face = grid.Face_Coord(axis, i);
					f(axis, face);
				}
			}
		}

		template<class IFFuncT>
		void Calc_Faces(IFFuncT f) {
			for (int axis = 0; axis < d; axis++) {
				const int dof = grid.Face_DoF(axis);
				thrust::counting_iterator<int> idxfirst(0);
				thrust::counting_iterator<int> idxlast = idxfirst + dof;
				thrust::transform(
					idxfirst,
					idxlast,
					face_data[axis]->begin(),
					[f, axis, this](const int idx) {
						return f(axis, grid.Face_Coord(axis, idx));
					}
				);
			}
		}
	};

	template<class T, int d> using FaceFieldDv = FaceField<T, d, DEVICE>;

}