<div align="center">


# flow2d-rs
Eulerian Fluid Simulation written in Rust.

Backward facing step

<img src="https://github.com/1n0r1/flow2d-rs/assets/80285371/49c0a938-7f22-43f3-a640-57d9d9d6be14" alt="backward facing step" width="1000"/>

Cylinder cross flow

<img src="https://github.com/1n0r1/flow2d-rs/assets/80285371/a4787740-4f53-4499-b82a-38897f698688" alt="cylinder cross flow" width="500"/>

Lid-driven cavity

<img src="https://github.com/1n0r1/flow2d-rs/assets/80285371/565b13a5-b87e-4e17-af37-9a7a0708b0c5" alt="lid-driven cavity" width="400"/>

</div>

## Features
- 2D viscous incompressible Newtonian fluid flow
- Solve Navier-Stokes equations using Euler's Method on staggered grid
- Solve Poisson equation with Successive Over-relaxation method
- Calculate and visualize pressure, speed and stream function
- Different types of boundary conditions:
  - Moving no-slip boundary condition
  - Free slip boundary condition
  - Inflow and outflow condition
- Planned features:
  - Standalone crate just for the simulation (without GUI)
  - Contour plot for stream function
  - Free boundary value simulation
  - HDF5 data export
  - Use other method to solve Poisson equation (possibly Multigrid)
  - Optimize to run on GPU
  - Energy/Heat flow simulation
  - Extension to 3D
- The theory and algorithm can be found in _Numerical simulation in fluid dynamics: a practical introduction_[[1]](#1)

## Quick start
```bash
  git clone https://github.com/1n0r1/flow2d-rs.git
  cd flow2d-rs
  cargo run
```
Refer to `./src/presets.rs` for setting up other simulations.

## Dependencies
- [Iced](https://github.com/iced-rs/iced) - for GUI
- [plotters](https://github.com/plotters-rs/plotters) - for exporting images


## References
<a id="1">[1]</a> Michael Griebel, Thomas Dornseifer, and Tilman Neunhoeffer. 1998. _Numerical simulation in fluid dynamics: a practical introduction_. SIAM
