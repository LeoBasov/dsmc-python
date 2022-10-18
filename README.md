# Octree based DSMC - Python implementation
[![LeoBasov](https://circleci.com/gh/LeoBasov/dsmc-python/tree/clean_up.svg?style=svg)](https://app.circleci.com/pipelines/github/LeoBasov/dsmc-python/tree/clean_up/)

# Table of Contents
1. Introduction
2. Requirements
3. Installation
4. Test Cases

# 1. Introduction
DSMC implementation using an octree as a variable mesh.
Implementation in done in python 3.10.

# 2. Requirements
- python 3.10.
- pip3

# 3. Installation
The module can be installed using pip3 from the repository root as

``
pip3 install .
``

# 4. Test Cases
All test cases were performed using Argon.
The gass properties are as follows:

| gas | $\sigma_T / m^2$ | $m / kg$   | 
|:---:|:----------------:|:----------:|
| Ar  | 3.631681e-19     | 6.6422e-26 |

## 4.1 Heat Bath
Simulation of temperature relaxation of Argon in closed domain.
The simulation domain is cube with a side length of $2 \cdot 10^{-3} m$.
The simulation properties are as follows


| $\Delta t / s$ | $w$    | $T / K$   | $n / m^{-3}$   | $u / (m/s)$ |
|:----------:|:------:|:-----:|:-----:|:-----------:|
| 1e-5       | 0.5e-8 | 300   | 1e+20 | 1000.0      |

where $\Delta t$ is the time step, $w$ is the particle weight, $T$  the temperature, $n$ the number density and $u$ the velocity in x direction.
The simulation results can be seen below.

![Heat Bath](./examples/heat_bath/heat_bath.png)

## 4.2 Hypersonic flow around cube
Hypersonic flow around a cuboid.
The inpflow conditions are as follows

| $T / K$ | $n / m^{-3}$ | $v_{x, z} / (m s^{-1})$ | $v_y / (m s^{-1})$ |
|:-------:|:------------:|:-----------------------:|:------------------:|
| 273.0   | 2.6e+19      | 0                       | -3043.0            |

![hypersonic_flow](./examples/hypersonic_flow/hypersonic_flow.png)

![hypersonic_flow_grid](./examples/hypersonic_flow/hypersonic_flow_grid.png)

![hypersonic_flow_grid_nrho](./examples/hypersonic_flow/hypersonic_flow_grid_nrho.png)

## 4.3 Shock Tube

ToDo
![shock Tube](./examples/shock_tube/shock_tube.png)
