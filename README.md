# N-Body Gravitational Simulation with Numerical Integrators

## Overview ------------------

This project is an interactive 2D N-body gravitational simulation written in C using Raylib.
It simulates a simplified solar system and compares different numerical time-integration schemes by analyzing energy conservation.

The main goal of this project is not visual realism, but to study numerical behavior of integrators used in classical mechanics.

## Motivation ------------------

In classical mechanics, most real systems do not have closed-form solutions.
Numerical integration is therefore essential, but different methods behave very differently over long time scales.

This project explores:

1. How numerical errors accumulate

2. Why higher-order methods are not always better

3. Why symplectic integrators are preferred for long-term Hamiltonian systems

## Features ------------------

1. Simulation of a 9-body gravitational system (Sun + 8 planets)

2. Interactive GUI to switch integrators

3. Real-time visualization of:

4. Orbits

5. Velocity vectors

6. Relative total energy error

7. Export of simulation data to CSV

8. Pause / play / reset functionality

## Numerical Integrators Implemented ------------------

1. Euler Method

    - First-order explicit integrator

    - Simple but highly unstable

    - Shows rapid energy drift

2. Symplectic Euler

    - First-order symplectic method

    - Much better long-term behavior than Euler

    - Energy error oscillates instead of drifting

3. RK2 (Midpoint Method)

    - Second-order Runge–Kutta

    - Better short-term accuracy

    - Still non-symplectic → long-term drift

4. RK4 (Classical Runge–Kutta)

    - Fourth-order accurate

    - Very good short-term precision

    - Does not conserve energy over long times

    - Floating-point round-off becomes visible at high magnification

5. Leapfrog (Velocity Verlet)

    - Second-order symplectic integrator

    - Best long-term energy behavior

    - Energy error remains bounded and oscillatory

## Physics Model ------------------

Newtonian gravity:

- F = GMm/|r|^3

- Softening term added to avoid singularities:

- r = sqrt(dx*dx + dy*dy + 1e-6)

- Initial conditions chosen to approximately satisfy:

- Circular orbits

- Zero total momentum

- Energy Diagnostics

The total energy is computed as:

𝐸 = 𝐾 + 𝑈

Where:

𝐾=∑1/2𝑚𝑣^2

𝑈=−∑𝑖<𝑗𝐺𝑚𝑖𝑚𝑗𝑟𝑖𝑗
​
The graph shows relative energy error:

𝐸(𝑡)−𝐸(0)/𝐸(0)​

- A large scaling factor is deliberately used to:

- Visualize tiny numerical errors

- Expose floating-point noise and long-term drift

## Key Observations

- Higher order ≠ better long-term stability

- RK4 performs extremely well short-term, but energy drifts slowly

- Symplectic methods preserve qualitative physics

- Leapfrog outperforms RK4 for long-time orbital simulations

- Numerical error behavior depends on both method and time step

## Controls

- Select integrator from menu

- Pause / Play simulation

- Toggle velocity vectors

- Reset simulation

- Rescale energy graph

- Change integrator at runtime

## Output

- Simulation data is saved as:

   Simulation_data.csv


- Each row contains:

  body_index, time, x, y, vx, vy, KE, PE, TE


This allows further analysis using Python / MATLAB / Excel.

## Dependencies

- C compiler (GCC recommended)

- Raylib

- raygui

## How to Run
- gcc main.c -o simulator -lraylib -lm
- ./simulator
(Adjust library paths if needed.)

## What This Project Demonstrates

- Practical implementation of numerical integrators

- Understanding of Hamiltonian dynamics

- Importance of symplectic methods

- Numerical error analysis

- Scientific visualization

## Future Improvements (Not Implemented)

- Adaptive time stepping

- Angular momentum diagnostics

- 3D extension

- Relativistic corrections

These were intentionally avoided to keep the focus on numerical behavior.

## Project-Specific Remarks

- The simulation uses a fixed time step

dt = 0.002


This value was chosen as a compromise between numerical stability and real-time performance.
Smaller time steps reduce energy drift but increase computational cost.

- For the purpose of visualizing numerical error, the energy error graph is deliberately scaled by a large factor.
This does not indicate unphysical behavior; it is done to make very small floating-point deviations visible on screen.

- Different integrators may therefore appear to produce “large” energy fluctuations on the graph, even when the actual relative error is extremely small.

- Forces are recomputed multiple times per step for RK2 and RK4, as required by their formulations.
This increases floating-point round-off effects, which become visible only because of the chosen visualization scale.

- All body masses remain constant throughout the simulation, ensuring that differences in behavior arise purely from the numerical integration method, not from changing physical parameters.

- The primary goal of this project is comparative numerical behavior, not high-precision ephemerides or astrophysical realism.

Author

Manas Kushwaha
IIT Kanpur
Mathematics & Scientific Computing

