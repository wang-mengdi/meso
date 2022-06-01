//////////////////////////////////////////////////////////////////////////
// Driver of two phase incompressible fluid solver on an Eulerian grid 
// Copyright (c) (2018-), Bo Zhu Shiying Xiong Zhecheng Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#define _USE_MATH_DEFINES
#include <math.h>
#include "SPX_Common.h"
#include "SPX_File.h"
#include "SPX_GeometryPrimitives.h"
#include "FluidEulerTwoPhaseClebsch.h"
#include "SPX_Driver.h"
#include "SPX_AuxFunc.h"
#include "RenderFunc.h"
#include "EulerInitializer.h"

using namespace AuxFunc;
template<int d> class FluidEulerTwoPhaseClebschDriver : public Driver
{Typedef_VectorDii(d);using Base=Driver;
public:
	FluidEulerTwoPhaseClebsch<d> fluid;
	real h_bar;
	real beta;
	real rho_A;
	real rho_L;
    bool use_bdry_v=true;
    bool use_velocity_field = false;

	FluidEulerTwoPhaseClebschDriver():fluid(){}

	real CFL() const
	{
		real epsilon=(real)1e-5;
		return (real)1*cfl*fluid.mac_grid.grid.dx/(fluid.Max_Abs(fluid.velocity)+epsilon);
	}

	virtual void Advance_One_Time_Step(const real dt,const real time)
	{
		fluid.Advance(dt);
	}

	virtual void Write_Output_Files(const int frame)
	{	
		Timer timer;timer.Reset();
		Base::Write_Output_Files(frame);
		if(frame==0){std::string file_name=frame_dir+"/grid";
			fluid.mac_grid.grid.Write_To_File_3d(file_name);
			std::cout<<"Write to file "<<file_name<<std::endl;}
		bool output_vel = true; ///v
		bool output_solid_phi = false;  ///S
		bool output_fluid_phi = true;  ///F
		bool output_psi = true;  ///L
		bool output_pressure = false;  ///P
		bool output_psi_N = true;  ///N
		bool output_psi_D = false;  ///D
		bool output_energy = true;
		
		////Write velocity
		if(output_vel==true){std::string file_name=frame_dir+"/velocity";
			fluid.velocity.Write_To_File_3d(file_name);}

		////Write psi
		if(output_psi){fluid.psi.Write_To_File_3d(frame_dir+"/psi_L");}

		////Write BC
		if(output_psi_D){RenderFunc::Write_Dirichlet_Boundary_Conditions(frame_dir+"/psi_D",fluid.mac_grid,fluid.bc);}
		if(output_psi_N){RenderFunc::Write_Neumann_Boundary_Conditions(frame_dir+"/psi_N",fluid.mac_grid,fluid.bc);}

		////Write fluid_phi
		if(output_fluid_phi) {
			fluid.levelset.phi.Write_To_File_3d(frame_dir + "/fluid_phi");
			//// levelset mesh
			RenderFunc::Write_LevelSet_Mesh(frame_dir + (d == 2 ? "/segment_mesh" : "/triangle_mesh"), fluid.levelset);
		}

		////Write energy
		if(output_energy){
			////Write energy
			real energy = fluid.Kinetic_Energy(fluid.velocity);
			real enstrophy = fluid.Enstrophy();
			real helicity = fluid.Helicity();
			std::cout << "Kinetic Energy for frame " << frame << ": " << energy << std::endl;
			std::cout << "Enstrophy for frame " << frame << ": " << enstrophy << std::endl;
			std::cout << "Helicity for frame " << frame << ": " << helicity << std::endl;
			std::string file_name = output_dir + "/0/energy.txt";
			std::string info = std::to_string(frame) + "\t" + std::to_string(energy) + "\t" + std::to_string(enstrophy) + "\t" + std::to_string(helicity) + "\n";
			if (frame == 0) { File::Write_Text_To_File(file_name, "Frame\tenergy\tenstrophy\thelicity\n"); }
			File::Append_Text_To_File(file_name, info);
		}

		timer.Elapse_And_Output_And_Reset("write to frame time");

		std::cout<<"Write to frame "<<frame<<std::endl;
	}
	
	virtual void Initialize()
	{
		int s=scale;real length=(real)1;VectorDi cell_counts=VectorDi::Ones();cell_counts[0]=s;
		cell_counts[1]=cell_counts[0]/2;if(d>2)cell_counts[2]=cell_counts[0]/2;
		switch(test){
		case 1: {	////hydrostatic/surface tension/gravity test
			length = 5;
			cell_counts = VectorDi::Ones() * s;
			fluid.Initialize(cell_counts, (real)length / cell_counts[0]);

			EulerInitializer<d> perimeter;
			perimeter.Set_Boundary_Width(0, 0, 0, 0, 0, 0);
			perimeter.Set_Boundary_Value(0, 0, 0, 0, 0, 0);
			perimeter.Set_Parameters(fluid.mac_grid);
			perimeter.Fill_Boundary_Condition(fluid.bc);

			fluid.projection.sigma = (real)1.e-2;

			////init levelset
			VectorD center = fluid.mac_grid.grid.Center();
			real r = fluid.mac_grid.grid.Length()[1] * (real).1;
			VectorD e_r = VectorD::Ones() * r; e_r[1] *= (real)2.;
			ImplicitShape<d> surface;
			//surface += std::make_shared<Plane<d>>(VectorD::Unit(1), center);						// hydrostatic
			surface += std::make_shared<Ellipsoid<d>>(center, e_r);							// surface tension
			//surface += std::make_shared<Sphere<d>>(center+0.5*center[1]*VectorD::Unit(1), r);		// gravity
			//fluid.use_body_force = true;
			fluid.levelset.Set_By_Shape(surface);
			fluid.levelset.Fast_Marching();

			//volume preserving
			//fluid.projection.use_vol_control = true;
			//fluid.projection.target_vol = fluid.levelset.Total_Volume();
			//fluid.projection.vol_control_ks = (real)1;

			fluid.initial_psi = std::bind(&FluidEulerTwoPhaseClebschDriver<d>::Initial_Zero, this, std::placeholders::_1);

			fluid.Update_Rho_Face();
			fluid.Update_Cell_Types();
		}break;
		case 10: { // horseshoe vortex
			 //./fluid_euler_two_phase_clebsch.exe -test 10 -lf 500 -s 64 -d 3 -cfl 1. -o horseshoe_64_clf_1_22
			length = 2;
			cell_counts = VectorDi::Ones() * scale;
			cell_counts[0] *= 2;
			if (d == 3) cell_counts[2] *= 2;
			fluid.Initialize(cell_counts, (real)length / cell_counts[0]);

			fluid.use_body_force = true;

			const auto mac_grid = fluid.mac_grid;

			EulerInitializer<d> perimeter;
			perimeter.Set_Boundary_Width(0, 0, 0, -1, 0, 0);
			perimeter.Set_Boundary_Value(0, 0, 0, 0, 0, 0);
			perimeter.Set_Parameters(mac_grid);
			perimeter.Fill_Boundary_Condition(fluid.bc);

			////Water bulk
			Plane<d> plane(VectorD::Unit(1), fluid.mac_grid.grid.Center()*4/3/* - VectorD::Unit(1) * (real)1*/);
			////Initialize phi
			fluid.levelset.Set_By_Geom(plane);

			////initialize water bulk psi
			VectorD bdry_v = VectorD::Unit(0) * (real)-0.;
			VectorD bdry_v_hbar = bdry_v / fluid.h_bar;
			fluid.initial_psi = std::bind(&FluidEulerTwoPhaseClebschDriver<d>::Initial_Vortex_Ring, this, std::placeholders::_1, bdry_v_hbar, length);

			fluid.Update_Cell_Types();
			fluid.Update_Rho_Face();
  		}break;
		default: {
			AuxFunc::Crash_With_Info("Please Select A Valid Cases");
		}break;
		}
		fluid.use_velocity_field = use_velocity_field;
		fluid.Initialize_Wave_Func();
		fluid.h_bar = h_bar;
		fluid.rho_A = rho_A;
		fluid.rho_L = rho_L;
		fluid.beta = 0.5;
		for (int i = 0; i < 5; i++) {
			fluid.Blend_Velocity();
			std::cout << "before projection" << std::endl; fluid.Divergence_Power();
			fluid.Enforce_Incompressibility((real)0.);
			fluid.Blend_Velocity();
			std::cout << "after projection" << std::endl; fluid.Divergence_Power();
		}
		fluid.beta = beta;
		//fluid.Extrapolation();
	}

	Vector4 Initial_Zero(const VectorD& pos)
	{
		Vector<C, 2> psi = Vel_To_Psi_C<d>(VectorD::Zero(), pos);
		return C2V(psi);
	}

	Vector4 Initial_Constant_Velocity(const VectorD& pos, const VectorD& bdry_v){
		Vector<C, 2> psi = Vel_To_Psi_C<d>(bdry_v, pos);
		return C2V(psi);
	}

	Vector4 Initial_Vortex_Ring(const VectorD& pos, const VectorD& bdry_v, const real length)
	{
		////Initial_Vortex_Ring
		Vector<C, 2> psi = Vel_To_Psi_C<d>(bdry_v, pos);
		VectorD c = V<d>((real)length / (real)3, (real)length / 3, (real)length / 2);
		real r = (real)0.2 * length;
		real rx = (pos[0] - c[0]) / r;
		real r2 = (pos - c).squaredNorm() / pow(r, 2);
		real fR = exp(-pow(r2 / (real)9, (real)4));
		C q(2. * rx * fR / (r2 + 1), (r2 + 1. - 2. * fR) / (r2 + 1));
		psi[0] *= q;
		return C2V(psi);
	}
};