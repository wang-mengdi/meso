# Taylor Vortex

The initial velocity of the field is created by two Taylor vortices.
The vorticity distribution of one vortex is given as:
$$
    \bm{w}(x,y)=\frac{U}{a}\left(2-\frac{r^2}{a^2}\right)\exp\left(\frac{1}{2}\left(1-\frac{r^2}{a^2}\right)\right)
$$
Here $U$ is the maximum tangential velocity, $a$ is the core size of the vortex and $r$ is the distance from the vortex to the given position. We use $U=1$ and $a=0.3$ here and two vortices are initialized at position $(-0.4,0)$ and $(0.4,0)$ in a $[-\pi,\pi]^2$ field.
We can take the sum of the for infulence from multiple vortices.

The correspondence of velocity and vorticity is as follows:
$$
  -\nabla\cdot\nabla\bm{u}=\nabla\times\bm{w}
$$
In order to solve for $\bm{u}$, we need to solve $\phi$ first.
$$
    -\Delta \bm{\phi}=\bm{w}
$$
Then we can solve for $\bm{u}$ using the curl operator.
$$
    \bm{u}=\nabla\times\bm{\phi}
$$

Ideally, we use periodic boundary condition. A good indicator of the performance of a fluid solver is that the enstrophy and the kinetic energy are conserved as the velocity field evolves.

reference:\
https://thesis.library.caltech.edu/2246/5/HolaMS.pdf