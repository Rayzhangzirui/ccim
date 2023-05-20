# Compact Coupling Interface Method (CCIM) for Elliptic Interface Problem

### Elliptic Interface Problem
The repository demonstrates the Compact Coupling Interface Method (CCIM) [1] for the following Elliptic Interface Problem in 3D. 

$$
\begin{cases}
-\nabla \cdot(\epsilon \nabla  u) + a u = f & \text{ in }\Omega \setminus \Gamma \\
\left[u\right] = \tau,\quad \left[\epsilon\nabla  u\cdot \mathbf{n} \right] = \sigma & \text{on }\Gamma \\
u = g & \text{on }\partial \Omega. 
\end{cases}
$$

Here $\Gamma$ is an interface that separates $\Omega \subset R^d$ into an inside region $\Omega^-$ and an outside region $\Omega^+$.  $\Gamma$ is represented implicitly by the zero level set of a function $\phi$ so that $\Gamma =\{ x| \phi(x) = 0\}$. $\mathbf{n}$ is the outward unit normal vector of the interface. g is the Dirichlet boundary condition. $\epsilon$, $f$, $a$: $\Omega \to R$ are given functions that might be discontinuous across $\Gamma$. The notation $[v]$ stands for jump of $v$ across the interface.  

### CCIM
CCIM obtains second order accurate solution and second order accurate gradient at the interface for complex interface in 3D

### Linear System Solver
CCIM results in an asymmetric sparse linear system. The code includes self-contained implementation of BICGSTAB (with ILU preconditioner) and Hypre[4] (currently without MPI).

### Level Set Method
The code also includes examples of a moving interface using the level set method
$$\phi_t + v_n |\nabla \phi| = 0$$
where the normal velocity of the interface is $v_n = [\nabla u \cdot \mathbb{n}]$, and $u$ is the solution of the elliptic interface problem.
We use the forward Euler method for time stepping, Godunov scheme for the Hamiltonian, WENO scheme, and the Fast marching Method to extend $v_n$.

### Other Finite Difference Method
Additionally, the repository also includes implementation of the Coupling Interface Method (CIM)and the Improved Coupling Interface Method (ICIM). 

### Usage



### Reference
[1] Zhang, Z., Cheng, L.-T., 2021. A Compact Coupling Interface Method with Accurate Gradient Approximation for Elliptic Interface Problems. https://doi.org/10.48550/arXiv.2110.12414

[2] Chern, I.-L., Shu, Y.-C., 2007. A coupling interface method for elliptic interface problems. Journal of Computational Physics 225, 2138–2174. https://doi.org/10.1016/j.jcp.2007.03.012

[3] Shu, Y.-C., Chern, I.-L., Chang, C.C., 2014. Accurate gradient approximation for complex interface problems in 3D by an improved coupling interface method. Journal of Computational Physics 275, 642–661. https://doi.org/10.1016/j.jcp.2014.07.017

[4] Falgout, R.D., Yang, U.M., 2002. hypre: A Library of High Performance Preconditioners, in: Sloot, P.M.A., Hoekstra, A.G., Tan, C.J.K., Dongarra, J.J. (Eds.), Computational Science — ICCS 2002, Lecture Notes in Computer Science. Springer, Berlin, Heidelberg, pp. 632–641. https://doi.org/10.1007/3-540-47789-6_66

