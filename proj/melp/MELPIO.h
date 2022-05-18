//////////////////////////////////////////////////////////////////////////
// Rendering functions
// Copyright (c) (2022-), Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Interpolation.h"
#include "Mesh.h"
#include "EulerParticles.h"

#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>

namespace Meso {
	template<int d>
	void Write_MELP_E(const EulerParticles<d>& e_particles, std::string file_name) {
		Typedef_VectorD(d); Typedef_MatrixD(d);
		const Array<VectorD>& pos = e_particles.xRef();
		const Array<VectorD>& vel = e_particles.uRef();
		const Array<MatrixD>& frame = e_particles.ERef();
		const Array<real>& nden = e_particles.ndenRef();
		// setup VTK
		vtkNew<vtkXMLUnstructuredGridWriter> writer;
		vtkNew<vtkUnstructuredGrid> unstructured_grid;

		vtkNew<vtkPoints> nodes;
		nodes->Allocate(pos.size());
		vtkNew<vtkDoubleArray> posArray;
		posArray->SetName("Position");
		posArray->SetNumberOfComponents(d);
		vtkNew<vtkDoubleArray> velArray;
		velArray->SetName("Velocity");
		velArray->SetNumberOfComponents(d);
		vtkNew<vtkDoubleArray> normArray;
		normArray->SetName("Normal");
		normArray->SetNumberOfComponents(d);
		vtkNew<vtkDoubleArray> ndenArray;
		ndenArray->SetName("NDen");
		ndenArray->SetNumberOfComponents(1);

		for (int i = 0; i < pos.size(); i++) {
			Vector3 pos3 = VectorFunc::V<3>(pos[i]);
			nodes->InsertNextPoint(pos3[0], pos3[1], pos3[2]);

			Vector3 vel3 = VectorFunc::V<3>(vel[i]);
			velArray->InsertNextTuple3(vel3[0], vel3[1], vel3[2]);

			Vector3 norm3 = VectorFunc::V<3>((VectorD)frame[i].col(d - 1));
			normArray->InsertNextTuple3(norm3[0], norm3[1], norm3[2]);

			ndenArray->InsertNextTuple1(nden[i]);
		}

		unstructured_grid->SetPoints(nodes);

		unstructured_grid->GetPointData()->AddArray(velArray);
		unstructured_grid->GetPointData()->AddArray(normArray);
		unstructured_grid->GetPointData()->AddArray(ndenArray);
		unstructured_grid->GetPointData()->SetActiveVectors("Velocity");

		writer->SetFileName(file_name.c_str());
		writer->SetInputData(unstructured_grid);
		writer->Write();
	}
}