# Adaptive-TFPS
This project compute the numerical solution of steady state radiative transfer equation discretized by DOM and Adaptive TFPS, which is capable of adaptively compress the angular domain using the information from optical parameters.

# Description of the code
new_trial.m: the script which computes the numerical solution of steady state radiative transfer equation discretized by DOM and Adaptive TFPS.

assemble.m: offline assembling of discrete operators and revealing low rank structures.

solution.m: online construction of right hand side vector, special solution and final solution.

HG.m: discretize the kernel function.

myfind.m: get the row index, column index and element value for all elements in a matrix.

qnwlege1.m, qnwlege1.m and qnwlege1.m: Discretie Ordinate Method for 1d, 2d and 3d case.

Figure/plot_figure.m: the script which plot the figures in our article.

# Description of the data files in the folder Figure
buffer_center_phiN6d1.mat: the scalar flux at cell centers computed by DOM and Adaptive TFPS with N=6 (84 discrete velocity directions), $\delta=10^{-1}$ for buffer zone case .

buffer_center_phiN6dinf.mat: the scalar flux at cell centers computed by DOM and Adaptive TFPS with N=6 (84 discrete velocity directions), $\delta=10^{-\infty}=0$ for buffer zone case .

buffer_vecsize_N2d1, buffer_vecsize_N2d3 and buffer_vecsize_N2d5: the number of Adaptive TFPS basis functions in each physical cell with N=2 (12 discrete velocity directions), $\delta=10^{-1}$, $10^{-3}$ or $10^{-5}$ for buffer zone case .

buffer_vecsize_N6d1, buffer_vecsize_N6d3 and buffer_vecsize_N6d5: the number of Adaptive TFPS basis functions in each physical cell with N=6 (84 discrete velocity directions), $\delta=10^{-1}$, $10^{-3}$ or $10^{-5}$ for buffer zone case .

lattice_center_phiN6d1.mat: the scalar flux at cell centers computed by DOM and Adaptive TFPS with N=6 (84 discrete velocity directions), $\delta=10^{-1}$ for lattice case .

lattice_center_phiN6dinf.mat: the scalar flux at cell centers computed by DOM and Adaptive TFPS with N=6 (84 discrete velocity directions), $\delta=10^{-\infty}=0$ for lattice case .

lattice_vecsize_N2d1, lattice_vecsize_N2d3 and lattice_vecsize_N2d5: the number of Adaptive TFPS basis functions in each physical cell with N=2 (12 discrete velocity directions), $\delta=10^{-1}$, $10^{-3}$ or $10^{-5}$ for lattice case.

lattice_vecsize_N6d1, lattice_vecsize_N6d3 and lattice_vecsize_N6d5: the number of Adaptive TFPS basis functions in each physical cell with N=6 (84 discrete velocity directions), $\delta=10^{-1}$, $10^{-3}$ or $10^{-5}$ for lattice case.

crho00, crho30, crho23: the value of $\mathop{\max}\limits_{\delta} \Vert E_{\delta,\mathfrak{i}}^{-1}\Vert_{2}$ for different choice of N (N=2, 3, 4 or 5) when $g_{C_{-}}=g_{C_{+}}=0$, $g_{C_{-}}=0.3, g_{C_{+}}=0$, or $g_{C_{-}}=0.2, g_{C_{+}}=-0.3$, respectively.

matrixA: the sparsity pattern for the linear operator discretized by DOM and TFPS with $I\times I=4\times 4$ grid and $M=3$.

# How to use the code to reproduce the results
Firstly run new_trial.m to get the data files in the folder Figure (or use the data files in the folder Figure directly). Then run plot_figure.m to get the figures in our article for RTE in xy geometry cases.
