<p align="center">
  <a href="https://github.com/philipplk/metis"><img alt="metis" src="logo.png" width="100%"></a>
  <p align="center">Computing Constrained Mechanical Systems.</p>
</p>

See [license](LICENSE) and [acknowledgements](#acknowledgements). Cite as:

```bibtex
@misc{kinon_metis_2021,
  author = {Kinon, Philipp Lothar},
  title = {metis: computing constrained dynamical systems (GitHub repository)},
  year = {2021},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/philipplk/metis}}
}
```

* * *

## Description

Metis was a mythical Titaness who became the goddess of wisdom, prudence and ingenuity. Ingenuity being the "practical, complex and implicit" kind differs from the other types of wisdom.

Likewise, this project targets the efficient and aesthetically pleasing numerical computation of dynamical systems being subject to holonomic constraints: particle systems as well as rigid body systems.

Metis is an object-oriented MATLAB code package - tested with the R2020b version.

A startscript loads the desired input-file. This input-file includes all necessary parameters for the given mechanical problem (geometry, loads, initial values,...), the chosen numerical integration scheme (time-step size, method, simulation time), the Postprocessing routine (plot quantities, animation, export) and the solution technique (max. iterations, tolerance). Metis creates all necessary objects and computes the approximate solution based on the given parameters. Eventually, one can choose to have an animation of the solution, some plots are created and the results are being exported.
The simulation is will be tracked by a log-file `./metis.log`.

#### Prerequisites

-   [MATLAB](https://www.mathworks.com/products/matlab.html)
-   [matlab2tikz][urlmatlab2tikz] (optional)


#### Theoretical Background

-   Initial Value Problems (IVP) of Constrained Dynamics leading to Differential-algebraic Equations (DAEs)
-   Numerical Integration (Direct Methods, Variational Integrators, Energy-Momentum Schemes)
-   Newton-Rhapson Method
-   Particle Systems
-   Rigid Body Dynamics with Director Formulation

[Back To The Top](#description)

* * *

## How To Start

1.  Clone this directory or download the .zip folder
2.  Get [matlab2tikz][urlmatlab2tikz] (optional)
3.  Open the MATLAB editor or run it with the shell script [metis.sh](metis.sh)
4.  Open [start_metis_single_analysis.m](start_metis_single_analysis.m)
5. Adjust `<input_file_name>` corresponding to a file from [/input](/input), for more info look at [README_input](/input/README_input.md)

```matlab
[simulation, system, integrator, solver] = Metis('input/<input_file_name>',1,1);
```
6. Adjust the path to the matlab2tikz directory in the chosen input-file
7.  Execute [start_metis_single_analysis.m](start_metis_single_analysis.m) for a first simulation
8.  Edit or change input file or create a new one in [/input](/input)
9.  For error analyses run [start_metis_error_analysis.m](start_metis_error_analysis.m) with corresponding input-file
10.  Have fun!

[Back To The Top](#description)

* * *

## References

Hamiltonian Dynamics:
- [Leimkuhler, Reich: Simulating Hamiltonian Dynamics, 2005](https://doi.org/10.1016/S0167-2789(99)00054-8)

Differential-Algebraic Equations:
- [Kunkel, Mehrmann: Differential-algebraic equations: analysis and numerical solution, 2006](https://doi.org/10.4171/017)

Numerical Integration:
- [Hairer, Lubich, Wanner: Geometric numerical integration: structure-preserving algorithms for ordinary differential equations, 2006](https://doi.org/10.1007/3-540-30666-8)
- [Lew, Mata: A Brief Introduction to Variational Integrators, 2016](https://doi.org/10.1007/978-3-319-31879-0_5)
- [Marsden, West: Discrete mechanics and variational integrators, 2001](https://doi.org/10.1017/S096249290100006X)
- [Leyendecker, Marsden, Ortiz: Variational integrators for constrained dynamical systems, 2008](https://doi.org/10.1002/zamm.200700173)
- [Gonzalez: Time integration and discrete Hamiltonian systems, 1996](https://doi.org/10.1007/BF02440162)

Rigid Body Dynamics with Directors:
- [Betsch, Steinmann: Constrained integration of rigid body dynamics, 2001](https://doi.org/10.1016/S0045-7825(01)00283-3)

[Back To The Top](#description)

* * *

## Author Info

-   GitHub - [philipplk](https://github.com/philipplk)
-   Website - [Philipp L. Kinon, Karlsruhe Institute of Technology (KIT)](https://www.ifm.kit.edu/english/14_5490.php)
-   ResearchGate - [Philipp L. Kinon](https://www.researchgate.net/profile/Philipp-Kinon)
-   LinkedIn - [Philipp Kinon](https://www.linkedin.com/in/philipp-kinon-9092781b5/)

[Back To The Top](#description)

* * *

## Acknowledgements

Metis was mainly initialised by Philipp L. Kinon during his master thesis at the [_Institute of Mechanics (IFM)_](https://www.ifm.kit.edu/english/index.php) at [_Karlsruhe Institute of Technology (KIT)_](https://www.kit.edu/english/), Germany. This project would not have been possible without the support of others:

**Academic Supervisors:**

-   [Prof. Dr.-Ing. Peter Betsch](https://www.ifm.kit.edu/english/14_4655.php)
-   [Simeon Schneider, M. Sc.](https://www.ifm.kit.edu/english/14_4890.php)
-   [Dr.-Ing. Mark Schiebl](https://www.ifm.kit.edu/english/14_4906.php)

**Coding support:**

[Julian Bauer M. Sc.](https://www.ifm.kit.edu/english/14_5166.php)
(GitHub: [JulianKarlBauer](https://github.com/JulianKarlBauer))

* * *

[Back To The Top](#description)

[urlmatlab2tikz]: https://github.com/matlab2tikz/matlab2tikz
