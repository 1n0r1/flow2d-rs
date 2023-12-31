<div align="center">


# flow2d-rs
Computational fluid dynamics library written in Rust


https://github.com/1n0r1/flow2d-rs/assets/80285371/36516766-5e4f-436e-8b6b-90dad9e28668


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
  cargo run --example gui
```
Refer to `./src/presets.rs` for setting up other simulations.

## Use this in your project
Add `flow2d_rs` as a dependency in your `Cargo.toml`

```toml
flow2d_rs = "0.1.0"
```

## Dependencies
- [Rayon](https://github.com/rayon-rs/rayon) - to parallelize computation


## References
<a id="1">[1]</a> Michael Griebel, Thomas Dornseifer, and Tilman Neunhoeffer. 1998. _Numerical simulation in fluid dynamics: a practical introduction_. SIAM
