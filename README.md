# Stabilizability Analysis with Polytopic Viability Kernel

## Abstract
This package is used to conduct failure detection for humanoid robot under multi-contact scenarios. After the specification of active contact status in **/user/hrp2/InitContact.txt**, a randmonized feasible robot's configuration will be generated and a robot world XML file will be produced as well for simulation. A default QP whole-body controller is chosen to stabilize the robot into a static configuration while 7 failure detection methods are implemented to monitor robot's state during robot's whole simulation process.

### Important files in (/src)
**Main.cpp**: This file contains main function. It initializes a feasible robot's state, generates a world XML file, and passes these prerequiste information for simulation.

**SimulationTest.cpp**: Inner simulation functions are included in this file. It simulates robot's performance step by step while calling functions to calculate robot's commanded configuration for stabilization and trajectories of failure detection result for 7 methods.

**StabilizingControllerContact.cpp**: This file contains functions to calculate robot's whole-body controller with a QP approach. Robot's acceleration, joint torques, magnitudes of contact forces, and violation of contact acceleration constraints are minimized in a weighted sum fashion. The optimized acceleration is then used to update robot's configuration.

**ConvexPolyhedron.cpp**: This file contains functions related to construction of contact polytope and support polygon, computation of fall indicators based on orbital energy, capture point, and Polytopic Viability Kernels.
* PIPGenerator(): This function genearates a number of planar inverted pendulums based on active contacts and outputs a fall indicator based on an assumption that robot's CoM moves in a fashion which stabilizes itself globally (HJB method). 
* RBGenerator(): This function calculates robot's fall indicator with a rigid-body assumption (Orbital Energy).
* CPCEGenerator(): This function calculates robot's fall indicator with a Capture Point Approach.
* ZeroStepCapturabilityGenerator(): This function evaluates robot's failure with an assumption that robot destabilizes itself in the same direction of its centroidal velocity at the instant of disturbance while satisfying frictional contact.
* ZMPGeneratorAnalysis(): Zero-moment point is used to for fall evaluation.

> Note that based on the different input (Contact Polytope or Support Polygon), **RBGenerator()** and **CPCEGenerator()** can output PVK-RB/CP, and OE/CP, respectively.

### Dependencies
Package compilation requres:
* [Klampt](https://github.com/krishauser/Klampt): An open-source, cross-platform software package for robot modeling, simulating, planning, optimization, and visualization.
* [SNOPT](https://github.com/snopt): A library for constrained optimization using Sequential Quadratic Programming method.
* [Gurobi](https://www.gurobi.com): A commercial optimization solver which is capable of solving LP, QP, MIQP, etc.
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page): An open-source header-only c++ library for linear algebra, and operations on vector, matrix and tensor.
* [Cddplus](https://github.com/cddlib/cddplus): A C++ implementation of the Double Description Method.
* [Cilantro](https://github.com/kzampog/cilantro): A lean C++ library for conduct convex hull algorithm.
