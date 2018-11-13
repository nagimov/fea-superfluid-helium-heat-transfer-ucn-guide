This code runs FEA analysis of helium-4 to helium-3 heat exchanger for UCN project [triumf.ca/UCN](http://www.triumf.ca/ucn).

Conventional FEA tools (e.g. Ansys) cannot be used without modifications of their solvers for this problem due to a unconventional character of heat transfer mechanisms:
* heat conductivity within superfluid helium-4 (a.k.a. helium-II) is described by [mutual friction theory](https://www.jstor.org/stable/100163) of superfluid helium
* heat transfer at liquid-solid interfaces at low temperatures is limited by [Kapitza resistance](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.41.48)

# System description

The system consists of superfluid liquid helium-4 contained in a cylindrical tube with adiabatic walls. One end of this tube is a subject of heat deposition from heavy radiation fields. Other end of the tube is surrounded by copper shell. Outer side of the copper shell is finned and immersed into colder liquid helium-3, forming a heat exchanger.

# Source code

Project consists of two parts: mesher (`mesher.py`) and solver (`solver.c`).

## Mesher

Mesher defines geometry, initial and boundary conditions of the system and generates axisymmetric mesh. In order to simplify solver code, mesher also creates list of FEA operations for the solver.

System geometry, mesh size and initial conditions are defined in `mesher.py`. List of FEA operations is generated within `mesher.py` using information obtained from generated mesh (e.g. geometry of finite elements) and external data (e.g. heat deposition).

Data is passed to the solver via csv files:

* `mesh_t.csv` contains initial temperature of the mesh
* `mesh_ops.csv` contains list of FEA operations describing heat transfer; operations include every possible combination of mesh pairs and necessary information for the finite heat transfer calculation (geometry and material codes)
* `heat_ops.csv` contains list of FEA operations describing heat deposition (beam heating); operations include list of mesh elements which are subject to heat deposition, as well as heat load per element

Sizes of mesh data array and FEA operations data arrays are saved to `params.h`.

## Solver

Solver iterates over the entire mesh using elementary mesh operations described in `mesh_ops.csv` and mesh temperature initialized in `mesh_t.csv`.

After completing the integration, solver saves final temperature distribution within the mesh to `mesh_t_post.csv`. Formats of `mesh_t_post.csv` and `mesh_t.csv` match, and `mesh_t_post.csv` can be re-used as initial conditions for the solver by renaming it to `mesh_t.csv`.

Thermal properties of helium-3, helium-4 and copper are defined in `thermal*.h`. There are three versions available:

* `thermal3.h` - cubic approximations
* `thermal2.h` - square and linear approximations
* `thermal1.h` - square, linear and constant approximations

Choice of `thermal*.h` affects the speed of integration, especially if thermal properties are set to be updated frequently (defined in `solver.c`).

Heat transfer equations between all possible combinations of materials are defined in `solver.c`. Integrator parameters are defined in `solver.c` as well.

In order to avoid use of heap and dynamic memory allocation (and those terrifying memory pointers), static data arrays are defined using data array sizes from `params.h`.

## Post-processing

Simple post-processing script can be used to display the temperature distribution after integration. In the provided example, `post.py` reads `mesh_t_post.csv` and displays it using matplotlib.

For monitoring purposes, `run_solver.py` script can be used to launch the compiled solver.
