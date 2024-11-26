# CoSTAR - Continuation of Solution Torus AppRoximations

© 2024 Engineering Dynamics Group, Institute of Mechanics, University of Kassel

CoSTAR (Continuation of Solution Torus AppRoximations) is a MATLAB toolbox developed by the Engineering Dynamics Group at the University of Kassel. It enables the efficient computation and analysis of stationary, periodic, and quasi-periodic solutions of dynamic systems. The toolbox is designed to calculate bifurcation diagrams using modular predictor-corrector methods that combine various approximation approaches, such as shooting and Fourier-Galerkin methods. This modular structure supports both research and application needs by facilitating method comparisons and problem-specific adaptations within a unified framework.

### **Key Features**
- **Solution Computation:** Includes equilibrium, periodic, and quasi-periodic solution methods. Currently, shooting methods and Fourier-Galerkin methods are implemented for periodic and quasi-periodic solutions.
- **Predictor-Corrector Framework:** A general curve-tracing algorithm that operates independently of the specific solution method used.
- **Stability Analysis:** Features under development include tools for stability evaluation via eigenvalue theory, monodromy matrices, Floquet multipliers, and Lyapunov exponents.

## **Getting Started**
To get started with CoSTAR, navigate to the [Tutorials folder](./Tutorials) in the repository. This folder contains a variety of tutorial scripts and live scripts that serve as an introduction to the toolbox. 

- All files in the `Tutorials` folder that begin with `Tutorial_` are recommended entry points.  
- Files with the `.m` extension provide standard MATLAB scripts, while `.mlx` files are MATLAB Live Scripts with an interactive format, making them ideal for exploration and visualization.

**Example Steps:**
1. Open MATLAB and navigate to the [Tutorials folder](./Tutorials).
2. Start with one of the `Tutorial_` files, such as `Tutorial_EQ.mlx`, for an introduction to the toolbox.
3. Follow the examples provided in the scripts to familiarize yourself with the key features and workflows.

### **Outlook**
Future developments will expand CoSTAR’s capabilities to track bifurcation points and create stability maps, further enhancing its utility for the analysis of complex dynamic systems.
