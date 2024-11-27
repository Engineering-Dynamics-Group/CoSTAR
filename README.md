# CoSTAR - Continuation of Solution Torus AppRoximations

© 2024 Engineering Dynamics Group, Institute of Mechanics, University of Kassel

CoSTAR (Continuation of Solution Torus AppRoximations) is a MATLAB toolbox developed by the Engineering Dynamics Group at the University of Kassel. It enables the efficient computation and analysis of stationary (equilibria, periodic and quasi-periodic) solutions of dynamic systems. The toolbox is designed to calculate bifurcation diagrams using modular predictor-corrector methods that combine various approximation approaches. This modular structure supports both research and application needs by facilitating method comparisons and problem-specific adaptations within a unified framework.

### **Key Features**
- **Solution Computation:** Equilibrium, periodic, and quasi-periodic (2 base frequencies) solutions can be computed based on the state-space representation of a dynamic system. Currently, finite-difference methods, Fourier-Galerkin methods (featuring error control) and (multiple) shooting methods are implemented for periodic and quasi-periodic solutions.
- **Continuation of Solution Branches:** Curve tracing is carried out by a general predictor-corrector algorithm that operates independently of the specific approximation method used. Various predictors, different step control methods and a live plot are available.
- **Stability Analysis:** Stability evaluation is available via eigenvalue theory (equilibrium solutions), Floquet multipliers obtained from the monodromy matrix (periodic solutions) and Lyapunov exponents (quasi-periodic solutions). Furthermore, common bifurcations points can be determined for equilibria and periodic solutions.
- **Help:** Tutorials, example scripts and a built-in help function provide comprehensive explanations and make it easy to get started with the toolbox.

For a more detailed list of all available features, please see the PDF [The MATLAB Toolbox CoSTAR](./The_MATLAB_Toolbox_CoSTAR.pdf).

### **Outlook**
The toolbox is continuously developed and future releases will expand CoSTAR’s capabilities. Possible upcoming improvements can be found in the PDF [The MATLAB Toolbox CoSTAR](./The_MATLAB_Toolbox_CoSTAR.pdf).

## **Getting Started**
To get started with CoSTAR, navigate to the [Tutorials folder](./Tutorials) in the repository. This folder contains a variety of tutorial and example scripts that serve as an introduction to the toolbox. 

- The Tutorials `Tutorials_EQ` (EQ: equilibria), `Tutorial_PS_...` (PS: periodic solutions) and `Tutorial_QPS_...` (QPS: quasi-periodic solutions) are recommended entry points.
- There is one tutorial for each solution type (`EQ`, `PS`, `QPS`) and approximation method (`FDM`: finite-difference method, `FGM`: Fourier-Galerkin method, `SHM`: (multiple) shooting method).
- Files with the `.m` extension provide standard MATLAB scripts, while `.mlx` files are MATLAB Live Scripts that feature MATLAB code alongside formatted text and equations.

**Example Steps:**
1. Open MATLAB and navigate to the [Tutorials folder](./Tutorials).
2. Start with one of the `Tutorial_...` files, such as `Tutorial_EQ.mlx`, for an introduction to the toolbox.
3. Follow the examples provided in the scripts to familiarize yourself with the key features and workflows.

If you want to have an overview of the code structure and the program flow of CoSTAR, please see the PDF [The MATLAB Toolbox CoSTAR](./The_MATLAB_Toolbox_CoSTAR.pdf).