#pragma once
namespace Meso
{
	enum CellType :char {
		INVALID = -1,
		INTERIOR,
		DIRICHLET,
		NEUMANN,
		BOUNDARY
	};
}