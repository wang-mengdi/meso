//////////////////////////////////////////////////////////////////////////
// Rendering functions
// Copyright (c) (2022-), Yitong Deng
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Interpolation.h"
#include "Mesh.h"
#include "SPHParticles.h"

//#include <vtkXMLStructuredGridWriter.h>
#include <vtkXMLPolyDataWriter.h>
//#include <vtkStructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkVertexGlyphFilter.h>

namespace Meso {

	template<int d>
	void Write_SPH(const SPHParticles<d>& particles, std::string file_name) {
		Typedef_VectorD(d); Typedef_MatrixD(d);
		const Array<VectorD>& pos = particles.xRef();
		const Array<VectorD>& vel = particles.uRef();
		const Array<real>& nden = particles.ndenRef();
		const Array<real>& rho = particles.rhoRef();
		const Array<int>& boundary = particles.BRef();
		//// setup VTK
		vtkNew<vtkXMLPolyDataWriter> writer;
		vtkNew<vtkPolyData> polyData;

		vtkNew<vtkPoints> nodes;
		nodes->Allocate(pos.size());
		vtkNew<vtkDoubleArray> posArray;
		posArray->SetName("Position");
		posArray->SetNumberOfComponents(3);
		vtkNew<vtkDoubleArray> velArray;
		velArray->SetName("Velocity");
		velArray->SetNumberOfComponents(3);
		vtkNew<vtkDoubleArray> ndenArray;
		ndenArray->SetName("NDen");
		ndenArray->SetNumberOfComponents(1);
		vtkNew<vtkDoubleArray> rhoArray;
		rhoArray->SetName("Rho");
		rhoArray->SetNumberOfComponents(1);
		vtkNew<vtkIntArray> boundaryArray;
		boundaryArray->SetName("Boundary");
		boundaryArray->SetNumberOfComponents(1);

		for (int i = 0; i < pos.size(); i++) {
			Vector3 pos3 = MathFunc::V<3>(pos[i]);
			nodes->InsertNextPoint(pos3[0], pos3[1], pos3[2]);

			Vector3 vel3 = MathFunc::V<3>(vel[i]);
			velArray->InsertNextTuple3(vel3[0], vel3[1], vel3[2]);

			ndenArray->InsertNextTuple1(nden[i]);
			boundaryArray->InsertNextTuple1(boundary[i]);
			rhoArray->InsertNextTuple1(rho[i]);
		}

		polyData->SetPoints(nodes);
		polyData->SetVerts(polyData->GetVerts());
		vtkNew<vtkVertexGlyphFilter> vertexGlyphFilter;

		polyData->GetPointData()->AddArray(velArray);
		polyData->GetPointData()->SetActiveVectors("velocity");

		polyData->GetPointData()->AddArray(ndenArray);
		polyData->GetPointData()->SetActiveScalars("nden");
		polyData->GetPointData()->AddArray(boundaryArray);
		polyData->GetPointData()->AddArray(rhoArray);

		writer->SetFileName(file_name.c_str());
		writer->SetInputData(polyData);
		writer->Write();
	}
}