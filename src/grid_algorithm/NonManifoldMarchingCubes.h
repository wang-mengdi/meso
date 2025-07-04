//////////////////////////////////////////////////////////////////////////
// Marching Cubes Algorithm
// Copyright (c) (2022-), Yunquan Gu
// This file is part of MESO, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#pragma once
#include "Grid.h"
#include "Field.h"
#include "Mesh.h"
#include <functional>

namespace Meso {
	using CellType = unsigned char; // support at most 256 types
	// Input:
	//		value_field: field that contains value from signed-distance-function
	//		label_field: field that indicates the label of vertex in filed
	// Output: 
	//		vertex_matrix: vertices
	//		element_matrix: edges
	// Note: Marching Square is a 2-d version of marching cubes. All inputs and outputs should admit a dimision of 2.
	//       This function is currently under construction.
	// Functionality:
	//		- Adapative: OK
	//		- Ambiguous resolve: N
	//		- Optimization: N
	//		- Support both Device and Host: N
	//		- Adapative interjunction: OK
	bool Point_Side(const Vector<real,2>& a, const Vector<real, 2>& b, const Vector<real, 2>& c) {
		return ((b.x() - a.x()) * (c.y() - a.y()) - (b.y() - a.y()) * (c.x() - a.x()))<0;
	}

	template<class T>
	void Non_Manifold_Marching_Square(
		VertexMatrix<T, 2>& vertex_matrix,
		ElementMatrix<2>& element_matrix,
		Field<std::function<CellType(Vector2)>,2>& find_label,
		const Field<CellType, 2, HOST>& label_field,
		const Field<T, 2, HOST>& value_field) 
	{
		Assert(label_field.grid.Is_Unpadded(), "Non_Manifold_Marching_Square field.grid {} padding not allowed", label_field.grid);
		Assert(value_field.grid.Is_Unpadded(), "Non_Manifold_Marching_Square field.grid {} padding not allowed", value_field.grid);

		// 0. Init
		const Grid<2>& grid = label_field.grid;							// const reference to access coordinates of grid
		const Vector2i cell_counts = grid.Counts() - Vector2i::Ones();	// compute cells counts by grid counts
		Array<Vector<T, 2>> vertices{}; Array<Vector2i> elements;		// vertices array is used to store generated vertices, similarly elements array is used to store lines

		// initilize vertex indice field for edge(x-axis and y-axis)
		// more specificly, the number stored in `vertex_index_on_edge` indicates the vertex index in `vertices`
		
		Field<int, 2> vertex_index_on_edge[2]; for (int i = 0; i < 2; i++) 
		{
			Array<Grid<2>> edge_grid(2); 
			for (int i = 0; i < 2; i++) {
				edge_grid[i] = Grid<2>(grid.Counts() - Vector2i::Unit(i), grid.dx);
				vertex_index_on_edge[i].Init(edge_grid[i], -1);
			}
		}

		//  1. Define basic function
		// `Hamming_Weight` is used to assist calculating set size in `Edge_Type`
		const auto Hamming_Weight = [](unsigned long long c) -> unsigned char
		{
			unsigned l = 0;
			for (; c != 0; l += 1) c &= c - 1;
			return l;
		};

		// `Edge_Type` takes arbitrary values of four labels and encode them into decimal value.
		// There are only 15 probabilities.
		const auto Edge_Type = [&Hamming_Weight](CellType v0, CellType v1, CellType v2, CellType v3) -> unsigned int
		{
			std::vector<CellType> v = { v0, v1, v2, v3 };
			unsigned char type = 0;
			for (size_t i = 0; i < 4; i++)
			{
				unsigned long long c = 0; // able to take 64 types of 
				for (size_t j = 0; j < i + 1; j++)
				{
					c |= (1ULL << v[j]);
					if (c & 1ULL << v[i])
						break;
				}
				type = type * 4 + Hamming_Weight(c);
			}
			return type - 85;
		};

		// `Gen_Vertex` take an edge(two point vector) as input and then generates vertex on it.
		// The generated vertex will be pushed into `vertices` and the `v_idx` will link the index with this edge.
		const auto Gen_Vertex = [&label_field, &value_field, &grid, &vertices](Field<int, 2>& v_idx, const Vector2i node_i, const Vector2i node_j)->void
		{
			const CellType& t1 = label_field(node_i); const CellType& t2 = label_field(node_j);
			const T& v1 = value_field(node_i); const T& v2 = value_field(node_j);

			if (t1 != t2)
			{
				T alpha = v1 * v2 < 0 ? (v1) / (v1 - v2) : (v1) / (v1 + v2);
				Vector2 pos = ((1 - alpha) * grid.Position(node_i) + alpha * grid.Position(node_j));
				vertices.push_back(pos.template cast<T>()); v_idx(node_i) = (int)vertices.size() - 1;
			};
		};

		// 2. Main part
		grid.Exec_Nodes(
			[&](const Vector2i cell_index) {
				// go through all cells and generate vertex on each four edges.
				if ((cell_index - cell_counts).maxCoeff() == 0) {
					find_label(cell_index) = [=](Vector2 pos) -> CellType {
						return label_field(cell_index);
					};
					return;
				}
				const CellType& left_down = label_field(cell_index);
				const CellType& right_down = label_field(cell_index + Vector2i::Unit(0));
				const CellType& right_up = label_field(cell_index + Vector2i::Ones());
				const CellType& left_up = label_field(cell_index + Vector2i::Unit(1));

				// Don't change the order!
				Gen_Vertex(vertex_index_on_edge[0], cell_index, cell_index + Vector2i::Unit(0)); // bottom edge
				Gen_Vertex(vertex_index_on_edge[1], cell_index + Vector2i::Unit(0), cell_index + Vector2i::Ones()); // right edge
				Gen_Vertex(vertex_index_on_edge[0], cell_index + Vector2i::Unit(1), cell_index + Vector2i::Ones()); // top edge
				Gen_Vertex(vertex_index_on_edge[1], cell_index, cell_index + Vector2i::Unit(1)); // left edge

				// in default, generate a vertex in the middle of cell
				// TODO: should be optimized later
				T inv_v1 = 1. / (-1e-8 + value_field(cell_index));
				T inv_v2 = 1. / (-1e-8 + value_field(cell_index + Vector2i::Unit(0)));
				T inv_v3 = 1. / (-1e-8 + value_field(cell_index + Vector2i::Unit(1)));
				T inv_v4 = 1. / (-1e-8 + value_field(cell_index + Vector2i::Ones()));

				Vector2 pos = ((grid.Position(cell_index) * inv_v1) \
					+ (grid.Position(cell_index + Vector2i::Unit(0)) * inv_v2) \
					+ (grid.Position(cell_index + Vector2i::Unit(1)) * inv_v3) \
					+ (grid.Position(cell_index + Vector2i::Ones())  * inv_v4)) / \
					(inv_v1 + inv_v2 + inv_v3 + inv_v4);
				vertices.push_back(pos.template cast<T>());

				const int center_vertex = (int)vertices.size() - 1;
				const int left_vertex = vertex_index_on_edge[1](cell_index);
				const int right_vertex = vertex_index_on_edge[1](cell_index + Vector2i::Unit(0));
				const int top_vertex = vertex_index_on_edge[0](cell_index + Vector2i::Unit(1));
				const int bottom_vertex = vertex_index_on_edge[0](cell_index);

				// Important: order: right-top -> left-top -> left-btm -> right-btm
				switch (Edge_Type(right_up, left_up, left_down, right_down))
				{
				// 0 segment
				case 0:
				{
					find_label(cell_index) = [=](Vector2 pos) -> CellType {
						return left_down;
					};
					break;
				}
				case 15:
				{
					find_label(cell_index) = [=](Vector2 pos) -> CellType {
						return left_down;
					};
					break;
				}
				// 1 segment
				case 1: 
				{
					elements.push_back(Vector2i(bottom_vertex, right_vertex));
					find_label(cell_index) = [=](Vector2 pos) -> CellType {
						return Point_Side(vertices[right_vertex], vertices[bottom_vertex], pos) ? left_down : right_down;
					};
					break;
				}
				case 4: 
				{
					elements.push_back(Vector2i(bottom_vertex, left_vertex));
					find_label(cell_index) = [=](Vector2 pos) -> CellType {
						return Point_Side(vertices[left_vertex], vertices[bottom_vertex], pos) ? left_down : right_down;
					};
					break;
				}
				case 5:
				{
					elements.push_back(Vector2i(left_vertex, right_vertex));
					find_label(cell_index) = [=](Vector2 pos) -> CellType {
						return Point_Side(vertices[left_vertex], vertices[right_vertex], pos) ? left_down : right_up;
					};
					break;
				}
				case 16:
				{
					elements.push_back(Vector2i(left_vertex, top_vertex));
					find_label(cell_index) = [=](Vector2 pos) -> CellType {
						return Point_Side(vertices[top_vertex], vertices[left_vertex], pos) ? left_up : right_down;
					};
					break;
				}
				case 20: {
					elements.push_back(Vector2i(top_vertex, bottom_vertex));
					find_label(cell_index) = [=](Vector2 pos) -> CellType {
						return Point_Side(vertices[top_vertex], vertices[bottom_vertex], pos) ? left_down : right_down;
					};
					break;
				}
				case 21: {
					elements.push_back(Vector2i(top_vertex, right_vertex));
					find_label(cell_index) = [=](Vector2 pos) -> CellType {
						return Point_Side(vertices[top_vertex], vertices[right_vertex], pos) ? left_down : right_up;
					};
					break;
				}
				// 2 segment
				case 6: 
				{
					elements.push_back(Vector2i(left_vertex, center_vertex));
					elements.push_back(Vector2i(center_vertex, right_vertex));
					elements.push_back(Vector2i(center_vertex, bottom_vertex));
					find_label(cell_index) = [=](Vector2 pos) -> CellType {
						return 0;
					};
					break;
				}
				case 17: // TODO: solve ambiguous
				{
					elements.push_back(Vector2i(left_vertex, top_vertex));
					elements.push_back(Vector2i(bottom_vertex, right_vertex));
					find_label(cell_index) = [=](Vector2 pos) -> CellType {
						return 0;
					};
					break;
				}
				case 18: 
				{
					elements.push_back(Vector2i(left_vertex, top_vertex));
					elements.push_back(Vector2i(bottom_vertex, right_vertex));
					find_label(cell_index) = [=](Vector2 pos) -> CellType {
						return 0;
					};
					break;
				}
				case 22: 
				{
					elements.push_back(Vector2i(top_vertex, center_vertex));
					elements.push_back(Vector2i(center_vertex, bottom_vertex));
					elements.push_back(Vector2i(center_vertex, right_vertex));
					find_label(cell_index) = [&](Vector2 pos) -> CellType {
						return 0;
					};
					break;
				}
				case 24: 
				{
					elements.push_back(Vector2i(top_vertex, center_vertex));
					elements.push_back(Vector2i(center_vertex, bottom_vertex));
					elements.push_back(Vector2i(center_vertex, left_vertex));
					find_label(cell_index) = [=](Vector2 pos) -> CellType {
						return 0;
					};
					break;
				}
				case 25: 
				{
					elements.push_back(Vector2i(left_vertex, bottom_vertex));
					elements.push_back(Vector2i(top_vertex, right_vertex));
					find_label(cell_index) = [=](Vector2 pos) -> CellType {
						return 0;
					};
					break;
				}
				case 26:
				{
					elements.push_back(Vector2i(left_vertex, center_vertex));
					elements.push_back(Vector2i(center_vertex, right_vertex));
					elements.push_back(Vector2i(top_vertex, center_vertex));
					find_label(cell_index) = [=](Vector2 pos) -> CellType {
						return 0;
					};
					break;
				}
				case 27:
				{
					elements.push_back(Vector2i(left_vertex, center_vertex));
					elements.push_back(Vector2i(center_vertex, right_vertex));
					elements.push_back(Vector2i(top_vertex, center_vertex));
					elements.push_back(Vector2i(center_vertex, bottom_vertex));
					find_label(cell_index) = [=](Vector2 pos) -> CellType {
						return 0;
					};
					break;
				}
				default:
				{
					//fmt::print("default: {} {} {} {} {}\n", right_up, left_up, left_down, right_down, Edge_Type(right_up, left_up, left_down, right_down));
					break;
				}
				}
			}
		);

		vertex_matrix = Eigen::Map<VertexMatrix<T, 2>>(vertices[0].data(), vertices.size(), 2);
		element_matrix = Eigen::Map<ElementMatrix<2>>(elements[0].data(), elements.size(), 2);
	}

	template<class T, int d, DataHolder side>
	void Non_Manifold_Marching_Cubes(VertexMatrix<T, d>& vertex_matrix, ElementMatrix<d>& element_matrix, Field<std::function<CellType(Vector<T,d>)>,d>& find_label, const Field<CellType, d, side>& label, const Field<T, d, side>& value) {
		Assert(label.grid.Is_Unpadded(), "marching cubes field.grid {} padding not allowed", label.grid);
		if constexpr (d == 2) {
			if constexpr (side == HOST) Non_Manifold_Marching_Square<T>(vertex_matrix, element_matrix, find_label, label, value);
			else Assert(false, "Marching_Squares not implemented for GPU");
		}
		else Assert(false, "Marching_Cubes not implemented for d={}", d);
	}

	template<class T, int d, DataHolder side>
	void Non_Manifold_Marching_Cubes(VertexMatrix<T, d>& vertex_matrix, ElementMatrix<d>& element_matrix, const Field<CellType, d, side>& label, const Field<T, d, side>& value) {
		Field<std::function<CellType(Vector<T,d>)>, d> find_label;
		find_label.Init(label.grid);

		Assert(label.grid.Is_Unpadded(), "marching cubes field.grid {} padding not allowed", label.grid);
		if constexpr (d == 2) {
			if constexpr (side == HOST) Non_Manifold_Marching_Square<T>(vertex_matrix, element_matrix, find_label, label, value);
			else Assert(false, "Marching_Squares not implemented for GPU");
		}
		else Assert(false, "Marching_Cubes not implemented for d={}", d);
	}
}