# 网格系统

三维按$4\times 4\times 4$切块，二维按$8\times 8$切块，这样每一个块内含有64个元素，适配CUDA的访存模式。

我们把网格分成CORNER和CENTER两种模式（Grid的模板参数grid_type）。直观上看，在CORNER模式下，节点（node）位于格子（cell）的角上，总共有`counts-VectorDi::Ones()`个格子；而在CENTER模式下，节点位于格子中心，总共有`counts`个格子。这里`counts`总是保存节点的数量。

另一个区别是，在CORNER模式下，第一个节点就位于domain_min。在CELL模式下，第一个格点是最小格子的中心，离domain_min在各个维度上均差0.5dx.

这种做法的意义是，对于一般的MAC网格，储存在cell中心的量构成了一个CENTER模式的Grid，而储存在面上的量，每一个轴向分别构成一个CORNER格式的Grid。如此便可把二者统一起来，即面上的数据亦为一个Grid，如此便不必分开处理cell数据和face数据。

domain_min这个参数仅在初始化的时候被使用，实际上会有一个`pos_min`变量，代表最小的节点位置，因为CORNER和CENTER模式的计算逻辑是一样的，这样不必每次都计算0.5dx，运算速度较快。

在CPX的设计中，Face_Index不同轴向所对应的编码方式不同。为简便起见，我们暂时忽略这种不同，一律按照axis==d-1的方式编码，即二维的axis=1，三维的axis=2，此时块内按照z-y-x的顺序编码（即，改z对应的index改变量最大，y次之，x最小）。若今后想要重新恢复这种编码方式的区分，可以通过给Grid类再加一个模板参数的方法实现。

# 网格系统的储存

我们保证，FaceField中每个轴向数据的储存顺序严格等同于grid使用Face_Grid()生成的CORNER模式grid的储存顺序。

Field和FaceField均按照“std::vector”的直觉设计，即所有的默认拷贝都是深拷贝，这样让用户不必操心未初始化指针的问题，可以放心复制。但实际上，其中的数据以std::shared_ptr形式存储，如此设计的原因是，可以临时地把FaceField在某个轴向上的数据取出，作为一个Field，执行Field专有的操作。
