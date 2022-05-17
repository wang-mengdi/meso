#include"levelset.h"

namespace Meso{
	/// Note !!!!
	/// Because: function template partial specialization is not allowed
	/// So I have to write code for DEVICE and HOST in one function
	/// which may make the code long.

template class LevelSet<2, PointIntpLinearPadding0,HOST>;
template class LevelSet<2, PointIntpLinearPadding0, DEVICE>;
template class LevelSet<3, PointIntpLinearPadding0, HOST>;
template class LevelSet<3, PointIntpLinearPadding0, DEVICE>;
//template class Grid<3, CELL>;
//template class Grid<3, NODE>;
}