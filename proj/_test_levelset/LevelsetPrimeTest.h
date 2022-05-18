//////////////////////////////////////////////////////////////////////////
// Test leveset on GPU and this is primary test
// Copyright (c) (2022-), Zhiqi Li
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "levelset.h"
#include "Timer.h"
using namespace Meso;
template<int d, class PointIntp,DataHolder side=DEVICE>
void Test_Circle(real dx, Vector<real, d> center = Vector<real, d>::Zero(), real radius = (real)1) {
	Typedef_VectorD(d);
	Sphere<d> geom(center, radius);
	
	//First test the init and assign
	VectorDi counts = (Vector<real, d>::Ones() * radius * 2.2 / dx).template cast<int>();
	Grid<d> grid(counts, dx, center - Vector<real, d>::Ones() * radius * 1.1);
	Timer timer;
	timer.Reset();
	LevelSet<d, PointIntp, HOST> levelset_host(grid);
	levelset_host.Set_By_Geom(geom);
	timer.Record("init with geom");
	LevelSet<d, PointIntp, side> levelset;
	levelset.Initialize(levelset_host);
	timer.Record("init with host");
	
	// begin to test the phi 
	Field<real, d, side> field_phi(grid);
	if constexpr(side==DEVICE) timer.Record("create device field_phi");
	else timer.Record("create host field_phi");
	levelset.All_Phi(field_phi);// levelset.*(&LevelSet<d, PointIntp, DEVICE>::Phi));
	timer.Record("calc field_phi");
	Array<real,HOST > field_phi2=*(field_phi.data);
	timer.Record("fetch field_phi");
	real maxError= 0;
	for (int i = 0; i < grid.DoF(); i++) {
		VectorDi cell = grid.Coord(i);
		VectorD pos = grid.Position(cell);
		real refPhi = geom.Phi(pos);
		real actPhi = field_phi2[i];
		if constexpr (d == 2) {
			Assert(fabs(refPhi - actPhi) < 1e-5, "Test_levelset_init failed: index {} failed,pos is {},{},ref phi is {},actual phi is {}", i, pos[0], pos[1], refPhi, actPhi);
		}
		else if constexpr (d == 3) {
			Assert(fabs(refPhi - actPhi) < 1e-5, "Test_levelset_init failed: index {} failed,pos is {},{},{},ref phi is {},actual phi is {}", i, pos[0], pos[1], pos[2], refPhi, actPhi);
		}
		if (fabs(refPhi - actPhi) > maxError) maxError = fabs(refPhi - actPhi);
	}
	Pass("Test_Prime_levelset_Phi (Circle) Passed, with error:{}", maxError);
	timer.Record("compare phi one by one");
	// then begin to test the normal
	Field<VectorD, d, side> field_normal(grid);
	if constexpr (side == DEVICE) timer.Record("create device field_normal");
	else timer.Record("create host field_normal");
	levelset.All_Normal(field_normal);// levelset.*(&LevelSet<d, PointIntp, DEVICE>::Phi));
	timer.Record("calc field_normal");
	Array<VectorD, HOST > field_normal2 = *(field_normal.data);
	timer.Record("fetch field_normal");
	maxError=0;
	for (int i = 0; i < grid.DoF(); i++) {
		//We not test the border
		VectorDi cell = grid.Coord(i);
		bool continue_flag = false;
		for (int j = 0; j < d; j++) {
			if (cell[j] == 0 || cell[j] == grid.counts[j] - 1) continue_flag = true;
		}
		if (continue_flag) continue;
		VectorD pos = grid.Position(cell);
		VectorD refNormal = geom.Normal(pos);
		VectorD actNormal = field_normal2[i];
		if constexpr (d == 2) {
			Assert((refNormal - actNormal).norm()/ actNormal.norm() < 0.05, "Test_levelset_Normal failed: index {} failed,pos is {},{},ref Normal is {},{},actual Normal is {},{}", i, pos[0], pos[1], refNormal[0], refNormal[1], actNormal[0], actNormal[1]);
		}
		else if constexpr (d == 3) {
			Assert((refNormal - actNormal).norm() / actNormal.norm() < 0.05, "Test_levelset_Normal failed: index {} failed,pos is {},{},{},ref Normal is {},{},{},actual Normal is {},{},{}", i, pos[0], pos[1], pos[2], refNormal[0], refNormal[1], refNormal[2], actNormal[0], actNormal[1], actNormal[2]);
		}
		if ((refNormal - actNormal).norm() / actNormal.norm() > maxError) maxError = (refNormal - actNormal).norm() / actNormal.norm();
	}
	Pass("Test_Prime_levelset_Normal (Circle) Passed, with error:{}",maxError);
	timer.Record("compare noraml one by one");
	// then begin to test the gridient
	Field<VectorD, d, side> field_gradient(grid);
	if constexpr (side == DEVICE) timer.Record("create device field_gradient");
	else timer.Record("create host field_gradient");

	levelset.All_Gradient(field_gradient);// levelset.*(&LevelSet<d, PointIntp, DEVICE>::Phi));
	timer.Record("calc field_gradient");
	Array<VectorD, HOST > field_gradient2 = *(field_gradient.data);
	timer.Record("fetch field_gradient");
	maxError = 0;
	for (int i = 0; i < grid.DoF(); i++) {
		//We not test the border
		VectorDi cell = grid.Coord(i);
		VectorD pos = grid.Position(cell);
		bool continue_flag = false;
		for (int j = 0; j < d; j++) {
			if (cell[j] == 0 || cell[j] == grid.counts[j] - 1) continue_flag = true;
			if ((pos - center).norm() < 2 * dx) continue_flag = true;
		}
		if (continue_flag) continue;
		
		VectorD refGradient = geom.Gradient(pos);
		VectorD actGradient = field_gradient2[i];
		if constexpr (d == 2) {
			Assert((refGradient - actGradient).norm() / actGradient.norm() < 0.1, "Test_levelset_Gradient failed: index {} failed,pos is {},{},ref Gradient is {},{},actual Gradient is {},{}", i, pos[0], pos[1], refGradient[0], refGradient[1], actGradient[0], actGradient[1]);
		}
		else if constexpr (d == 3) {
			Assert((refGradient- actGradient).norm() / actGradient.norm() < 0.1, "Test_levelset_Gradient failed: index {} failed,pos is {},{},{},ref Gradient is {},{},{},actual Gradient is {},{},{}", i, pos[0], pos[1], pos[2], refGradient[0], refGradient[1], refGradient[2], actGradient[0], actGradient[1], actGradient[2]);
		}
		if ((refGradient - actGradient).norm() / actGradient.norm() > maxError) maxError = (refGradient - actGradient).norm() / actGradient.norm();
	}
	Pass("Test_Prime_levelset_Gradient (Circle) Passed, with error:{}", maxError);
	timer.Record("compare gradient one by one");
	// then begin to test the curvature
	Field<real, d, side> field_curvature(grid);
	if constexpr (side == DEVICE) timer.Record("create device field_curvature");
	else timer.Record("create host field_curvature");
	levelset.All_Curvature(field_curvature);// levelset.*(&LevelSet<d, PointIntp, DEVICE>::Phi));
	timer.Record("calc field_curvature");
	Array<real, HOST > field_curvature2 = *(field_curvature.data);
	timer.Record("fetch field_curvature");
	maxError = 0;
	for (int i = 0; i < grid.DoF(); i++) {
		//We not test the border
		VectorDi cell = grid.Coord(i);
		VectorD pos = grid.Position(cell);
		bool continue_flag = false;
		for (int j = 0; j < d; j++) {
			if (cell[j] <=1 || cell[j] >= grid.counts[j] - 2) continue_flag = true;
			if ((pos - center).norm() < 4 * dx) continue_flag = true;
		}
		if (continue_flag) continue;

		real refCurvature = geom.Curvature(pos);
		real actCurvature = field_curvature2[i];
		if constexpr (d == 2) {
			Assert(fabs(refCurvature - actCurvature)/ fabs(actCurvature)< 0.1, "Test_levelset_Curvature failed: index {} failed,pos is {},{},ref Curvature is {},actual Curvature is {}", i, pos[0], pos[1], refCurvature, actCurvature);
		}
		else if constexpr (d == 3) {
			Assert(fabs(refCurvature - actCurvature)/ fabs(actCurvature) < 0.1, "Test_levelset_Curvature failed: index {} failed,pos is {},{},{},ref Curvature is {},actual Curvature is {}", i, pos[0], pos[1], pos[2], refCurvature, actCurvature);
		}
		if (fabs(refCurvature - actCurvature) / fabs(actCurvature) > maxError) maxError = fabs(refCurvature - actCurvature) / fabs(actCurvature);
	}
	Pass("Test_Prime_levelset_Curvature (Circle) Passed, with error:{}", maxError);
	timer.Record("compare Curvature one by one");
	timer.Output_Profile();
	if constexpr (side==DEVICE)
		Pass("Test_Prime_levelset_GPU (Circle) Passed!");
	else
		Pass("Test_Prime_levelset_CPU (Circle) Passed!");
}