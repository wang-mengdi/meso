//////////////////////////////////////////////////////////////////////////
// Rendering functions
// Copyright (c) (2022-), Mengdi Wang
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Mesh.h"
#include "AuxFunc.h"
#include <vtkNew.h>
#include <vtkTriangle.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>

using namespace Meso;
namespace VTKFunc {
	template<int d,class T>
	void Output_VTU(const std::shared_ptr< SurfaceMesh<d> >& mesh,const Array<T>& velocity, std::string file_name) {
		Typedef_VectorD(d);

		// setup VTK
		vtkNew<vtkXMLUnstructuredGridWriter> writer;
		vtkNew<vtkUnstructuredGrid> unstructured_grid;

		vtkNew<vtkPoints> nodes;
		nodes->Allocate(velocity.size());
		vtkNew<vtkDoubleArray> velArray;
		velArray->SetName("Velocity");
		velArray->SetNumberOfComponents(d);

		vtkNew<vtkCellArray> cellArray;

		for (int i = 0; i < mesh->Vertices().size(); i++) {
			Vector3 pos3 = VectorFunc::V<3>(mesh->Vertices()[i]);
			nodes->InsertNextPoint(pos3[0], pos3[1], pos3[2]);

			Vector3 vel3 = VectorFunc::V<3>(velocity[i]);
			velArray->InsertNextTuple3(vel3[0], vel3[1], vel3[2]);
		}

		for (int i = 0; i < mesh->Elements().size(); i++) {
			vtkNew<vtkTriangle> triangle;
			triangle->GetPointIds()->SetId(0, mesh->Elements()[i][0]);
			triangle->GetPointIds()->SetId(1, mesh->Elements()[i][1]);
			triangle->GetPointIds()->SetId(2, mesh->Elements()[i][2]);
			cellArray->InsertNextCell(triangle);
		}
		
		unstructured_grid->SetPoints(nodes);
		unstructured_grid->SetCells(VTK_TRIANGLE, cellArray);

		unstructured_grid->GetPointData()->AddArray(velArray);
		unstructured_grid->GetPointData()->SetActiveVectors("velocity");

		writer->SetFileName(file_name.c_str());
		writer->SetInputData(unstructured_grid);
		writer->Write();
	}
}
