# RoboDojo.jl
[![CI](https://github.com/thowell/RoboDojo.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/thowell/RoboDojo.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/thowell/RoboDojo.jl/branch/main/graph/badge.svg?token=E0EBU6ZJRL)](https://codecov.io/gh/thowell/RoboDojo.jl)

A differentiable simulator for robotic systems. This repository includes models for planar hopper, biped, and quadruped systems. 

Systems are represented in minimal coordinates and contact impulses are computed at each time step by solving a [nonlinear complementarity problem (NCP)](https://en.wikipedia.org/wiki/Nonlinear_complementarity_problem) using a bespoke [path-following solver](src/solver/interior_point.jl). Gradients through the dynamics are efficiently computed via the [implicit-function theorem](https://en.wikipedia.org/wiki/Implicit_function_theorem) applied to the optimality conditions of the NCP.

For more details, see our paper: [Fast Contact-Implicit Model-Predictive Control](https://arxiv.org/abs/2107.05616) and the related [repository](https://github.com/thowell/ContactImplicitMPC.jl).

If you are interested in using this package for your research or there are features you would like implemented, please reach out to [Taylor](thowell@stanford.edu) or [Simon](simon-lc@stanford.edu).

## Installation
``` 
git clone https://github.com/thowell/RoboDojo.jl
```

## Systems
### Hopper 

### Biped 

### Quadruped 

## Coming Soon 
- [ ] documentation 
- [ ] system schematics
- [ ] 3D quadruped model 
- [ ] Atlas model
- [ ] tutorial

