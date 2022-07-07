#pragma once
#include "Grid.h"
#include "Field.h"

namespace Meso {
	template<class T, int d>  Field<T, d> Convolution_Filter(real r, const Field<T, d>& f_in)
	{
		Typedef_VectorD(d);
		Field<T, d> f_out = f_in;
		const Grid<d>& grid = f_out.grid;
		int n = (int)(std::floor(r / f_in.grid.dx));
		Grid<d> sub_grid(VectorDi::Ones() * (n * 2 + 1), grid.dx);
		VectorDi sub_offset = VectorDi::Ones() * n;

		f_out.Calc_Nodes(
			[&](const VectorDi node)->real {
				real w = (real)0; real sum = (real)0;
				sub_grid.Iterate_Nodes(
					[&](const VectorDi sub_node)->void{
						const VectorDi& nb_node = sub_node + node - sub_offset;
						if (f_out.grid.Valid(nb_node)) {
							real d0 = (grid.Position(node) - grid.Position(nb_node)).norm();
							real w0 = std::max(r - d0, (real)0);
							sum += w0 * f_in(nb_node);
							w += w0;
						}
					}
				);
				if (w != (real)0) { return sum / w; }
				else { return f_in(node); }
			}
		);

		return f_out;
	}
}