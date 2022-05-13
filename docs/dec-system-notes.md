# 外微分系统的数学表示

### 微分形式

*本文主要内容摘抄自：《流体力学的量子化》，熊诗颖*

切空间上的$(0,r)$型反对称张量称为流形上的$r$形式：
$$
\begin{equation}
    t =\sum_{i_1<\cdots<i_r}t^{i_1,\cdots,i_r}\textrm{d}x^{i_1}\wedge\cdots \wedge \textrm{d}x^{i_r},
\end{equation}
$$


对于二维流形, 假设其上具有局部坐标系$x^1,x^2$, 则0形式 $\alpha$, 1形式 $\beta$, 和2形式 $\gamma$可以写成
$$
\begin{equation}
\begin{dcases}
    \alpha = \alpha,\\
    \beta = \beta_1 \textrm{d} x^1 + \beta_2 \textrm{d} x^2,\\
    \gamma = \gamma_{1,2} \textrm{d} x^1 \wedge \textrm{d} x^2,\\
\end{dcases}
\end{equation}
$$
在二维欧氏空间中，可以简单理解$x^1=x,x^2=y$.

对于三维流形, 假设其上具有局部坐标系$x^1,x^2,x^3$, 则0形式$\alpha$, 1形式$\beta$, 2形式$\gamma$, 和 3形式 $\delta$可以写成
$$
\begin{equation}
\begin{dcases}
    \alpha = \alpha,\\
    \beta = \beta_1 \textrm{d} x^1 + \beta_2 \textrm{d} x^2 + \beta_3 \textrm{d} x^3 ,\\
    \gamma = \gamma_{1,2} \textrm{d} x^1 \wedge \textrm{d} x^2+ \gamma_{2,3} \textrm{d} x^2 \wedge \textrm{d} x^3+ \gamma_{3,1} \textrm{d} x^3 \wedge \textrm{d} x^1,\\
    \delta = \delta_{1,2,3}\textrm{d} x^1 \wedge \textrm{d} x^2\wedge \textrm{d} x^3.
\end{dcases}
\end{equation}
$$
在三维欧氏空间中，可以简单理解$x^1=x,x^2=y,x^3=z$.

$r$形式是反对称张量，这意味着两个特性：

$$
\begin{equation}
\begin{dcases}
    \textrm{d}x^i\wedge\textrm{d}x^i=0,\\
    \textrm{d}x^i\wedge\textrm{d}x^j=-\textrm{d}x^j\wedge\textrm{d}x^i.
\end{dcases}
\end{equation}
$$

例如，在二维欧氏空间中：

$$
\begin{equation}
\begin{dcases}
    \textrm{d}x\wedge\textrm{d}x=\textrm{d}y\wedge\textrm{d}y=0,\\
    \textrm{d}x\wedge\textrm{d}y=-\textrm{d}y\wedge\textrm{d}x.
\end{dcases}
\end{equation}
$$

可以证明，$r$形式在进行坐标变换（即更换不同的基）时，其值按照相应的雅可比行列式进行放缩，这意味着$r$形式是一个可以积分的量，因此我们也将$r$形式叫做微分形式。

例如，在欧氏空间当中，我们可以简单地理解为，$\iint_S\textrm{d}x\wedge\textrm{d}y=\iint_S\textrm{d}x\textrm{d}y$，$\iiint_V\textrm{d}x\wedge\textrm{d}y\wedge\textrm{d}z=\iiint_V\textrm{d}x\textrm{d}y\textrm{d}z$.

将积分写成微分形式的版本, 除了坐标变换变得自然合理外, 还能统一斯托克斯, 高斯, 格林等公式, 即假设$\mathcal{M}$是一个$n$维流形, $u$是$n-1$阶形式, 则有
$$
\begin{equation}
    \int_{\mathcal{M}} \textrm{d} u = \int_{\partial \mathcal{M}} u
\end{equation}
$$
如果我们将速度场看成一形式, 则涡量场$\omega = \textrm{d} u$为二形式, 上式的直接含义就是涡量在一个曲面的积分 (即涡通量) 等于速度在其边界的积分 (即速度环量).

### 外微分算子$\textrm{d}$

现定义外微分算子$\textrm{d}$.它会把一个$r$形式转换为一个$r+1$形式，例如，对上一节中的$r$形式$t$有
$$
\begin{equation}
   \textrm{d}t =  \sum_{i_1<\cdots<i_r}\frac{\partial t_{i_1,\cdots,i_r}}{\partial x^j}\textrm{d}x^{j}\wedge\textrm{d}x^{i_1}\wedge\cdots \wedge \textrm{d}x^{i_r}.
\end{equation}
$$


举一个例子，对三维欧氏空间中的一个$2$形式$a\textrm{d}x\wedge\textrm{d}y$，我们可以得到

$$
\textrm{d}\left(a\textrm{d}x\wedge\textrm{d}y\right)=\frac{\partial a}{\partial x}\textrm{d}x\wedge\textrm{d}x\wedge\textrm{d}y+\frac{\partial a}{\partial y}\textrm{d}y\wedge\textrm{d}x\wedge\textrm{d}y+\frac{\partial a}{\partial z}\textrm{d}z\wedge\textrm{d}x\wedge\textrm{d}y.
$$

右端前两项均为零（因$\textrm{d}x$或$\textrm{d}y$出现了两次），而又注意到$\textrm{d}z\wedge\textrm{d}x\wedge\textrm{d}y=-\textrm{d}x\wedge\textrm{d}z\wedge\textrm{d}y=\textrm{d}x\wedge\textrm{d}y\wedge\textrm{d}z$，故结果为

$$
\frac{\partial a}{\partial z}\textrm{d}x\wedge\textrm{d}y\wedge\textrm{d}z.
$$

当然，外微分算子也可以应用于不同阶形式的组合。对上一节中二维流形的$0\sim2$阶形式求外微分，可得：
$$
\begin{equation}
\begin{dcases}
    \textrm{d}\alpha =\frac{\partial \alpha}{\partial x^1}\textrm{d} x^1 + \frac{\partial \alpha}{\partial x^2}\textrm{d} x^2,\\
    \textrm{d}\beta =\textrm{d} \beta_1 \wedge \textrm{d} x^1 +\textrm{d} \beta_2 \wedge  \textrm{d} x^2=\left(\frac{\partial \beta_2}{\partial x^1}-\frac{\partial \beta_1}{\partial x^2}\right)\textrm{d} x^1\wedge \textrm{d} x^2,\\
    \textrm{d}\gamma = \textrm{d}\gamma_{1,2} \wedge(\textrm{d} x^1 \wedge \textrm{d} x^2)=0.\\
\end{dcases}
\end{equation}
$$

同理，对三维流形的$0\sim 3$阶形式求外微分：

$$
\begin{equation}
\begin{dcases}
    \textrm{d}\alpha = \frac{\partial \alpha}{\partial x^1}\textrm{d} x^1 + \frac{\partial \alpha}{\partial x^2}\textrm{d} x^2+ \frac{\partial \alpha}{\partial x^3}\textrm{d} x^3,\\
    \textrm{d}\beta =\left(\frac{\partial \beta_2}{\partial x^1}-\frac{\partial \beta_1}{\partial x^2}\right)\textrm{d} x^1\wedge \textrm{d} x^2+\left(\frac{\partial \beta_3}{\partial x^2}-\frac{\partial \beta_2}{\partial x^3}\right)\textrm{d} x^2\wedge \textrm{d} x^3+\left(\frac{\partial \beta_1}{\partial x^3}-\frac{\partial \beta_3}{\partial x^1}\right)\textrm{d} x^3\wedge \textrm{d} x^1,\\
    \textrm{d}\gamma = \left(\frac{\partial \gamma_{2,3}}{\partial x^1}+\frac{\partial \gamma_{3,1}}{\partial x^2}+\frac{\partial \gamma_{1,2}}{\partial x^3} \right)\textrm{d} x^1 \wedge \textrm{d} x^2\wedge\textrm{d} x^3,\\
    \textrm{d}\delta =0
\end{dcases}
\end{equation}
$$

可以看出$\textrm{d}$对3维流形中0形式$\alpha$, 1形式$\beta$, 2形式$\gamma$的作用类似于向量微积分中的梯度$\nabla$, 旋度$\nabla\times$, 和散度$\nabla\cdot$. 


### 霍奇星算子$*$

$n$维流形中，霍奇星算子$*$可以把一个$r$形式转换为一个$n-r$形式。特别地，$n$维欧氏空间的度规张量是$n$阶单位矩阵，那么此时，霍奇星算子的放缩系数为$1$，可以理解为一个简单的转换。例如，在三维欧氏空间中，$*1=\textrm{d}x\wedge\textrm{d}y\wedge\textrm{d}z$，$*\textrm{d}x=\textrm{d}y\wedge\textrm{d}z$，$*(\textrm{d}y\wedge\textrm{d}z)=\textrm{d}x$，$*(\textrm{d}x\wedge\textrm{d}y\wedge\textrm{d}z)=1$.

霍奇星算子是使用外微分术语描述向量微积分的必备工具。例如，从上一节中我们得知，对$2$形式求$\textrm{d}$算子可得到一个$3$形式的散度，但事实上，在向量微积分当中，散度应当作用于一个向量，也就是$1$形式，而散度的结果是一个标量，也就是$0$形式。这就需要使用$*$算子进行转换。具体到这里，就是先用$*$把$1$形式转换为$2$形式，再做$\mathrm{d}$求出一个$3$形式的散度，再用$*$把$3$形式转回$0$形式。

还有一点需要注意，流形中的向量可以表示为基的线性组合，以三维欧氏空间的标准正交基为例，即$\bm{u}=u\bm{e}_x+v\bm{e}_y+w\bm{e}_z$，但我们需要的是$1$形式。欧氏空间当中比较简单，可以直接认为$\bm{e}_x=\mathrm{d}x,\bm{e}_y=\mathrm{d}y,\bm{e}_z=\mathrm{d}z$，那么$\bm{u}=u\mathrm{d}x+v\mathrm{d}y+w\mathrm{d}z$.

因此，三维欧氏空间中，向量微积分的各个算子实际对应如下：

$$
\begin{equation}
\begin{dcases}
\nabla \phi=\textrm{d}\phi,\quad \text{0-form}\rightarrow\text{1-form},\\
\nabla\cdot\bm{u}=*\mathrm{d}*\bm{u},\quad\text{1-form}\rightarrow\text{0-form},\\
\nabla\times\bm{u}=*\mathrm{d}\bm{u},\quad\text{1-form}\rightarrow\text{1-form},\\
\nabla^2\phi=*\mathrm{d}*\mathrm{d}\phi,\quad\text{0-form}\rightarrow\text{0-form}.
\end{dcases}
\end{equation}
$$

其中$\phi$为标量场，即$0$形式。可以把$\phi$和$\bm{u}=u\mathrm{d}x+v\mathrm{d}y+w\mathrm{d}z$代入验算此四式成立。

#MAC网格上的外微分系统

MAC网格中有两种存数据的位置，一种是格子中心，一种是格子的面心。我们可以发现，格子是一个三维的体积元，因此可以对应$3$形式（因为$3$形式可以在三维体积元中做积分），而类似地，面对应$2$形式。

至于$1$和$2$形式，则需要用到对偶的概念。使用类似图论中对偶图的概念，在对偶图中，MAC网格的一个格子对应一个节点，而MAC网格两个格子之间的面对应这两个节点之间的连线。节点对应$0$形式，而连线对应$1$形式。

因此，在格子中心储存的数据可以对应$0$或$3$形式，而在面上储存的数据可以对应$1$或$2$形式，它们之间互为对偶，这种对偶实际上就是霍奇星符号。从物理的角度，我们可以理解为：格子中心的数据可以是一个标量（$0$形式），也可以是一个标量在格子所占体积里的积分（$3$形式）；面上的数据可以是一个垂直于面的矢量（$1$形式），也可以是该矢量在面上积分所得的通量（$2$形式）。

二维空间的情况类似，格子储存$0$和$2$形式，面上储存$1$形式。

在MAC网格的程序当中，我们实现的$\mathrm{d}$算子名为`Exterior_Derivative()`函数。它有两种重载方式：

```c++
//input: C is a 0-form on cell (2d or 3d)
//output: F is a 1-form on face (3d) or 1-form on edge (2d)
template<class T, int d> void Exterior_Derivative(FaceFieldDv<T, d>& F, const FieldDv<T, d>& C);
//input: F is a 2-form on face(3d) or 1-form on edge (2d)
//output: C is a 3-form on cell (3d) or 2-form on cell (2d)
template<class T, int d> void Exterior_Derivative(FieldDv<T, d>& C, const FaceFieldDv<T, d>& F);
```

第一种是从$0$形式到$1$形式（也就是从格子到面），第二种是从$2$形式到$3$形式（也就是从面到格子）。它们分别对应不同形式下的离散$\mathrm{d}$算子，这里离散的意思是不考虑网格尺寸$\Delta x$，或认为$\Delta x=1$。

我们并没有实现霍奇星算子，因为在MAC网格的离散化当中，对偶操作仅仅代表改变对数据的解释，数据本身的值并未改变。

### 泊松系统的外微分实现

接下来我们以泊松系统为例，研究如何使用外微分方法实现向量微分算子。

代码中的`MaskedPoissonMapping`类实现了泊松系统。具体来说，它是一个从格子到格子的映射，若某个格子的mask（即成员变量`fixed`）为true，则对该格子作单位映射（即值不变）；否则对改格子作负的泊松映射$-\nabla\cdot(\beta\nabla p)$，其中$\beta$就是成员变量`vol`.

单位映射是为了处理边界条件。在一般的流体系统当中，固体和空气格不应该进入泊松系统求解，但我们不希望破坏网格的完整性，因此把它们作一单位映射。负泊松映射的原因是，一个正定的线性系统对数值求解较为有利，但泊松算子$\nabla^2$的离散形式是一个负定矩阵，加个符号就可以把它变成正定系统。

去除mask部分，`MaskedPoissonMapping::Apply`函数的核心操作可以总结如下：

```c++
//d(p) ----- 1-form
Exterior_Derivative(temp_face, temp_cell);

//d(p) *. vol ----- 1-form
temp_face *= vol;

//Hodge star is identity here
//*d(p) *. vol ----- 2-form

//temp_cell = div(temp_face)
//d*d(p) *. vol ----- 3-form
Exterior_Derivative(temp_cell, temp_face);
temp_cell *= -1;

//Hodge star is identity here
//-*d*d(p) *. vol ----- 0-form
```
其中，输入是$p$，而泊松算子对应$*\mathrm{d}*\mathrm{d}p$.对函数`Exterior_Derivative()`的两次调用就是其中的两个$\mathrm{d}$算子，而霍奇星是隐含的。

### 投影算法的外微分实现

另一个例子是MAC网格不可压流体的投影算法。`FluidEuler`类实现了这一功能，相关代码如下：

```c++
//projection
//vel_div=div(velocity)
Exterior_Derivative(vel_div, velocity);

int iter; real res;
MGPCG.Solve(pressure.Data(), vel_div.Data(), iter, res);
Info("Solve poisson with {} iters and residual {}", iter, res);

//velocity+=grad(p)
Exterior_Derivative(temp_velocity, pressure);
temp_velocity *= poisson.vol;

velocity += temp_velocity;
psi_N.Apply(velocity);
```

可以看到，求解泊松系统之前，从速度场$\bm{u}$求散度$\nabla\cdot\bm{u}$就是用了一个`Exterior_Derivative()`，也就是$\mathrm{d}$算子来完成。而求解泊松系统之后，计算压强$p$的梯度，亦采用`Exterior_Derivative()`函数完成。这种采用外微分语言的做法，和手动计算相比，极大降低了程序的复杂度。