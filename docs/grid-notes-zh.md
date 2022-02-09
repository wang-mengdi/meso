# 网格系统

三维按$4\times 4\times 4$切块，二维按$8\times 8$切块，这样每一个块内含有64个元素，适配CUDA的访存模式。

网格分成NODE和CELL模式（Grid的模板参数grid_type）。这二者的区别是它们和domain_min的关系。在NODE模式下，第一个格点就位于domain_min。在CELL模式下，第一个格点是最小格子的中心，离domain_min在各个维度上均差0.5dx.

domain_min这个参数仅在初始化的时候被使用，实际上会有一个`pos_min`变量，代表最小格点的位置，因为NODE和CELL模式的计算逻辑是一样的，这样不必每次都计算0.5dx，运算速度较快。
