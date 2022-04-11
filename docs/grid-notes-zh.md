# 网格系统

三维按$4\times 4\times 4$切块，二维按$8\times 8$切块，这样每一个块内含有64个元素，适配CUDA的访存模式。

我们把网格分成CORNER和CENTER两种模式（Grid的模板参数grid_type）。直观上看，在CORNER模式下，节点（node）位于格子（cell）的角上，总共有`counts-VectorDi::Ones()`个格子；而在CENTER模式下，节点位于格子中心，总共有`counts`个格子。这里`counts`总是保存节点的数量。

另一个区别是，在CORNER模式下，第一个节点就位于domain_min。在CELL模式下，第一个格点是最小格子的中心，离domain_min在各个维度上均差0.5dx.


domain_min这个参数仅在初始化的时候被使用，实际上会有一个`pos_min`变量，代表最小的节点位置，因为CORNER和CENTER模式的计算逻辑是一样的，这样不必每次都计算0.5dx，运算速度较快。

在CPX的设计中，Face_Index不同轴向所对应的编码方式不同。为简便起见，我们暂时忽略这种不同，一律按照axis==d-1的方式编码，即二维的axis=1，三维的axis=2，此时块内按照z-y-x的顺序编码（即，改z对应的index改变量最大，y次之，x最小）。若今后想要重新恢复这种编码方式的区分，可以通过给Grid类再加一个模板参数的方法实现。
