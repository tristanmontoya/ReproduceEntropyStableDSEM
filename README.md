# Efficient Entropy-Stable Discontinuous Spectral-Element Methods Using Tensor-Product Summation-by-Parts Operators on Triangles and Tetrahedra

This repository contains the Julia code to reproduce the results in the following manuscript:

T. Montoya and D. W. Zingg, "Efficient Entropy-Stable Discontinuous Spectral-Element Methods Using Tensor-Product Summation-by-Parts Operators on Triangles and Tetrahedra", In preparation, 2023.

Please cite the above manuscript if you use this repository or the underlying spectral-element framework [StableSpectralElements.jl](https://github.com/tristanmontoya/StableSpectralElements.jl) in your research. Any questions regarding the content of the manuscript or technical issues with the code should be directed to the corresponding author at [tristan.montoya@mail.utoronto.ca](mailto:tristan.montoya@mail.utoronto.ca).

## Abstract
We present a new class of efficient and robust discontinuous spectral-element methods of arbitrary order for nonlinear hyperbolic systems of conservation laws on curved triangular and tetrahedral unstructured grids. Such discretizations employ a recently introduced family of sparse tensor-product summation-by-parts (SBP) operators in collapsed coordinates within an entropy-conservative modal formulation, which is rendered entropy stable when a dissipative numerical flux is used at element interfaces. The proposed algorithms exploit the structure of such SBP operators alongside that of the Proriol-Koornwinder-Dubiner polynomial basis used to represent the numerical solution on the reference triangle or tetrahedron, and a weight-adjusted approximation is employed in order to efficiently invert the local mass matrix for curvilinear elements. Using such techniques, we obtain an improvement in time complexity from $\mathcal{O}(p^{2d})$ to $\mathcal{O}(p^{d+1})$ relative to existing entropy-stable formulations using multidimensional SBP operators not possessing such a tensor-product structure, where $p$ is the polynomial degree of the approximation and $d$ is the number of spatial dimensions. The number of required entropy-conservative two-point flux evaluations between pairs of quadrature nodes is accordingly reduced by a factor ranging from 1.52 at $p=2$ to 4.52 at $p=10$ for triangles, and from 1.88 at $p=2$ to 10.99 at $p=10$ for tetrahedra. Through numerical experiments involving smooth solutions to the compressible Euler equations on isoparametric triangular and tetrahedral grids, the proposed methods using tensor-product SBP operators are shown to exhibit similar levels of accuracy for a given mesh and polynomial degree to those using multidimensional operators based on symmetric quadrature rules, with both approaches achieving order $p+1$ convergence with respect to the element size in the presence of interface dissipation as well as exponential convergence with respect to the polynomial degree. Furthermore, both operator families are shown to give rise to entropy-stable schemes which exhibit excellent robustness for test problems characteristic of under-resolved turbulence simulations. Such results suggest that the algorithmic advantages resulting from the use of tensor-product operators are obtained without compromising accuracy or robustness, enabling the efficient extension of the benefits of entropy stability to higher polynomial degrees than previously considered for triangular and tetrahedral elements.

## Installation
First, make you have [Julia](https://julialang.org/downloads/) installed you haven't already done so. The tests in this paper were run on v1.9.3, and we recommend using that version or a later one. Then, assuming that you are using Linux or macOS and have git installed, follow the steps below.

1. Clone this repository by entering the command `git clone https://github.com/tristanmontoya/ReproduceEntropyStableDSEM.git` in the terminal.

2. Within the top-level `ReproduceEntropyStableDSEM` directory, use the command `julia --project=.` to open the Julia REPL and activate the project within the current directory. 

3. Install all dependencies by entering `using Pkg; Pkg.instantiate()` in the REPL. This will automatically set up the latest version of [StableSpectralElements.jl](https://github.com/tristanmontoya/StableSpectralElements.jl) for you to use within this project.

## Reproducibility instructions
Here, we describe how to generate the results using the provided drivers, and how to produce the results in the manuscript using the provided Jupyter notebooks. Note that the tests run significantly faster with multithreading enabled (for example, add `--threads 8` to the `julia` command if you want to use eight threads) as this allows for local element-based operations to be performed simultaneously. If using multiple Julia threads, it is [usually best to set the number of BLAS threads to 1](https://carstenbauer.github.io/ThreadPinning.jl/dev/explanations/blas/) (for example, using the `OPENBLAS_NUM_THREADS` environment variable). For all tests, enter the top-level directory and run `julia --project=.` and then load the driver package by entering `using ReproduceEntropyStableDSEM` in the REPL.

### Accuracy tests
To run the $h$-refinement studies on triangles, we set the variable `p` to the polynomial degree of the scheme (for example, the figures in the paper use `p=4` and `p=5`) and set the variable `results_dir` to the desired path for storing the results (for example, to store the results within this repository for use with the provided notebooks, we can take `results_dir="./results/"`). We can then run the following code snippets:
```julia
run_driver(EulerDriver(p, l=p, p_quad=35, C_t=0.5,
    element_type="Tri", 
    scheme="ModalTensor", 
    form="FluxDifferencingForm",
    numerical_flux="EntropyConservativeNumericalFlux", 
    ode_algorithm="DP8", 
    path=string(results_dir, "euler_h_refine/"), 
    M0=2, L=2.0, n_periods=1, mesh_perturb=1/16, n_grids=7, 
    load_from_file=true, overwrite=false, test_type=2))
```

```julia
run_driver(EulerDriver(p, l=p, p_quad=35, C_t=0.5,
    element_type="Tri", 
    scheme="ModalTensor", 
    form="FluxDifferencingForm",
    numerical_flux="LaxFriedrichsNumericalFlux", 
    ode_algorithm="DP8", 
    path=string(results_dir, "euler_h_refine/"), 
    M0=2, L=2.0, n_periods=1, mesh_perturb=1/16, n_grids=7, 
    load_from_file=true, overwrite=false, test_type=2))
```
```julia
run_driver(EulerDriver(p, l=p, p_quad=35, C_t=0.5,
    element_type="Tri", 
    scheme="ModalMulti", 
    form="FluxDifferencingForm",
    numerical_flux="EntropyConservativeNumericalFlux", 
    ode_algorithm="DP8", 
    path=string(results_dir, "euler_h_refine/"), 
    M0=2, L=2.0, n_periods=1, mesh_perturb=1/16, n_grids=7, 
    load_from_file=true, overwrite=false, test_type=2))
```

```julia
run_driver(EulerDriver(p, l=p, p_quad=35, C_t=0.5,
    element_type="Tri", 
    scheme="ModalMulti", 
    form="FluxDifferencingForm",
    numerical_flux="LaxFriedrichsNumericalFlux", 
    ode_algorithm="DP8", 
    path=string(results_dir, "euler_h_refine/"), 
    M0=2, L=2.0, n_periods=1, mesh_perturb=1/16, n_grids=7, 
    load_from_file=true, overwrite=false, test_type=2))
```
Likewise, we can run the $p$-refinement studies as follows:

```julia
run_driver(EulerPRefinementDriver(2,15, p_quad=35, C_t=0.5,
    element_type="Tri", 
    scheme="ModalTensor", 
    form="FluxDifferencingForm", 
    numerical_flux="LaxFriedrichsNumericalFlux", 
    ode_algorithm="DP8", 
    path=string(results_dir, "euler_p_refine/"), 
    M0 = 4, T=2.0, mesh_perturb = 1/16, 
    load_from_file=true, overwrite=false, test_type=1))
```
```julia
run_driver(EulerPRefinementDriver(2,15, p_quad=35, C_t=0.5,
    element_type="Tri", 
    scheme="ModalTensor", 
    form="FluxDifferencingForm", 
    numerical_flux="LaxFriedrichsNumericalFlux", 
    ode_algorithm="DP8", 
    path=string(results_dir, "euler_p_refine/"),
    M0 = 4, T=2.0, mesh_perturb = 1/16, 
    load_from_file=true, overwrite=false, test_type=1))
```
```julia
run_driver(EulerPRefinementDriver(2,15, p_quad=35, C_t=0.5,
    element_type="Tri", 
    scheme="ModalMulti", 
    form="FluxDifferencingForm", 
    numerical_flux="EntropyConservativeNumericalFlux", 
    ode_algorithm="DP8", 
    path=string(results_dir, "euler_p_refine/"),
    M0 = 4, T=2.0, mesh_perturb = 1/16, 
    load_from_file=true, overwrite=false, test_type=1))
```

```julia
run_driver(EulerPRefinementDriver(2,15, p_quad=35, C_t=0.5,
    element_type="Tri", 
    scheme="ModalMulti", 
    form="FluxDifferencingForm", 
    numerical_flux="EntropyConservativeNumericalFlux", 
    ode_algorithm="DP8",
    path=string(results_dir, "euler_p_refine/"),
    M0 = 4, T=2.0, mesh_perturb = 1/16, 
    load_from_file=true, overwrite=false, test_type=1))
```
The same studies can be repeated for tetrahedra case by replacing `element_type="Tri"` with `element_type="Tet"`, where the results in the manuscript use a smaller time step, replacing `C_t=0.5` with `C_t=0.15` in the tetrahedral case. For the $p$-refinement studies on tetrahedra, the multidimensional operators considered in our work are only up to $p=10$, in which case the corresponding function calls are:
```julia
run_driver(EulerPRefinementDriver(2,10, p_quad=35, C_t=0.15,
    element_type="Tet", 
    scheme="ModalMulti", 
    form="FluxDifferencingForm", 
    numerical_flux="EntropyConservativeNumericalFlux", 
    ode_algorithm="DP8", 
    path=string(results_dir, "euler_p_refine/"),
    M0 = 4, T=2.0, mesh_perturb = 1/16, 
    load_from_file=true, overwrite=false, test_type=1))
```

```julia
run_driver(EulerPRefinementDriver(2,10, p_quad=35, C_t=0.15,
    element_type="Tet", 
    scheme="ModalMulti", 
    form="FluxDifferencingForm", 
    numerical_flux="EntropyConservativeNumericalFlux", 
    ode_algorithm="DP8", 
    path=string(results_dir, "euler_p_refine/"),
    M0 = 4, T=2.0, mesh_perturb = 1/16, 
    load_from_file=true, overwrite=false, test_type=1))
```

### Robustness tests
The Kelvin-Helmholtz tests on triangles are run for a range of polynomial degrees from 4 to 9 using the following function calls for entropy-stable schemes using Lax-Friedrichs interface dissipation:
```julia
run_driver(EulerPRefinementDriver(4,9, C_t=0.025, 
    element_type="Tri",
    scheme="ModalTensor",
    form="FluxDifferencingForm",
    numerical_flux="LaxFriedrichsNumericalFlux",
    ode_algorithm="DP8", 
    path=string(results_dir, "kelvin_helmholtz/"),
    M0=16, T=15.0, mesh_perturb = 1/16, 
    load_from_file=true, overwrite=false, test_type=2))
```
```julia
run_driver(EulerPRefinementDriver(4,9, C_t=0.025, 
    element_type="Tri",
    scheme="ModalMulti",
    form="FluxDifferencingForm",
    numerical_flux="LaxFriedrichsNumericalFlux",
    ode_algorithm="DP8", 
    path=string(results_dir, "kelvin_helmholtz/"),
    M0=16, T=15.0, mesh_perturb = 1/16, 
    load_from_file=true, overwrite=false, test_type=2))
```
Likewise, the inviscid Taylor-Green vortex problem can be run with $\mathrm{Ma} = 0.1$ as follows: 
```julia
run_driver(EulerPRefinementDriver(4,9, C_t=0.005, 
    element_type="Tet",
    scheme="ModalTensor",
    form="FluxDifferencingForm",
    numerical_flux="LaxFriedrichsNumericalFlux",
    ode_algorithm="DP8",
    path=string(results_dir, "inviscid_tgv_ma01/"), 
    mach_number=0.1,
    M0=4, T=14.0, mesh_perturb = 1/16, 
    load_from_file=true, overwrite=false, test_type=2))
```
```julia
run_driver(EulerPRefinementDriver(4,9, C_t=0.005, 
    element_type="Tet",
    scheme="ModalMulti",
    form="FluxDifferencingForm",
    numerical_flux="LaxFriedrichsNumericalFlux",
    ode_algorithm="DP8", 
    path=string(results_dir, "inviscid_tgv_ma01/"),
    mach_number=0.1,
    M0=4, T=14.0, mesh_perturb = 1/16, 
    load_from_file=true, overwrite=false, test_type=2))
```

The Mach 0.7 cases can be run by substituting `mach_number=0.1` for `mach_number=0.7` in the above. Additionally, the cases can be run with an entropy-conservative numerical flux by substituting `numerical_flux="LaxFriedrichsNumericalFlux"` for  `numerical_flux="EntropyConservativeNumericalFlux"` as in the accuracy tests.

### Reproducing manuscript figures

The following table lists the Jupyter notebooks within the `notebooks` subdirectory which must be run in order to generate each figure in the manuscript.

|Figure| Description | Notebook(s)|
|---|---|---|
1 | 2D collapsed coordinate transformation  | [plot_collapsed_2d.ipynb](https://github.com/tristanmontoya/ReproduceEntropyStableDSEM/tree/main/notebooks/plot_collapsed_2d.ipynb) 
2 | 3D collapsed coordinate transformation  | [plot_collapsed_3d.ipynb](https://github.com/tristanmontoya/ReproduceEntropyStableDSEM/tree/main/notebooks/plot_collapsed_3d.ipynb) 
3 | Quadrature nodes for SBP operators | [plot_nodes_tri.ipynb](https://github.com/tristanmontoya/ReproduceEntropyStableDSEM/tree/main/notebooks/plot_nodes_tri.ipynb) </br>[plot_nodes_tet.ipynb](https://github.com/tristanmontoya/ReproduceEntropyStableDSEM/tree/main/notebooks/plot_nodes_tet.ipynb) 
4 | Number of entropy-conservative flux evaluations | [two_point_fluxes_tri.ipynb](https://github.com/tristanmontoya/ReproduceEntropyStableDSEM/tree/main/notebooks/two_point_fluxes_tri.ipynb) </br>[two_point_fluxes_tet.ipynb](https://github.com/tristanmontoya/ReproduceEntropyStableDSEM/tree/main/notebooks/two_point_fluxes_tet.ipynb) 
5 | Number of floating-point operations | [matrix_ops_tri.ipynb](https://github.com/tristanmontoya/ReproduceEntropyStableDSEM/tree/main/notebooks/matrix_ops_tri.ipynb) </br>[matrix_ops_tet.ipynb](https://github.com/tristanmontoya/ReproduceEntropyStableDSEM/tree/main/notebooks/matrix_ops_tet.ipynb) 
6 | Example meshes and mapping nodes | [plot_mesh_tri.ipynb](https://github.com/tristanmontoya/ReproduceEntropyStableDSEM/tree/main/notebooks/plot_mesh_tri.ipynb) </br>[plot_mesh_tet.ipynb](https://github.com/tristanmontoya/ReproduceEntropyStableDSEM/tree/main/notebooks/plot_mesh_tet.ipynb) 
7 | Convergence with respect to $h$ and $p$ | [accuracy_tri_h_refine.ipynb](https://github.com/tristanmontoya/ReproduceEntropyStableDSEM/tree/main/notebooks/accuracy_tri_h_refine.ipynb) </br>[accuracy_tri_p_refine.ipynb](https://github.com/tristanmontoya/ReproduceEntropyStableDSEM/tree/main/notebooks/accuracy_tri_p_refine.ipynb) </br>[accuracy_tet_h_refine.ipynb](https://github.com/tristanmontoya/ReproduceEntropyStableDSEM/tree/main/notebooks/accuracy_tet_h_refine.ipynb) </br>[accuracy_tet_p_refine.ipynb](https://github.com/tristanmontoya/ReproduceEntropyStableDSEM/tree/main/notebooks/accuracy_tet_p_refine.ipynb)
8 | Robustness tests (Kelvin-Helmholtz instability and Taylor-Green vortex) | [robustness_tri_khi.ipynb](https://github.com/tristanmontoya/ReproduceEntropyStableDSEM/tree/main/notebooks/robustness_tri_khi.ipynb) </br>[robustness_tet_tgv_ma01.ipynb](https://github.com/tristanmontoya/ReproduceEntropyStableDSEM/tree/main/notebooks/robustness_tet_tgv_ma01.ipynb) </br>[robustness_tet_tgv_ma07.ipynb](https://github.com/tristanmontoya/ReproduceEntropyStableDSEM/tree/main/notebooks/robustness_tet_tgv_ma07.ipynb)

 For Figures 7 and 8, the error norms and entropy histories required to generate the figures are provided in the `results` subdirectory, although the raw simulation datasets are not included in this repository due to their large file sizes. Such datasets can be generated by running the code as described above or obtained from the authors upon request. The figures generated using such notebooks are identical to those appearing in the manuscript, and are provided in the `plots` directory.

## License

This software is released under the [GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.en.html).