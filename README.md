# Recursive-polynomial-CAM-design
This repository includes all the necessary routines to run the CAM design with the method described in the paper in the articles [1,2].

The codes can be used to compute CAMs by employing any combination of the following features:

Multiple finite-burn dynamics and multi-impulsive dynamics
Earth Orbits dynamics, with the inclusion of second order harmonics
Circular Restricted Three-Body Problem in Cislunar
Multiple short-term encounters between the primary and more than one secondary
Grid filtering technique for selecting the most impactful firing times
Fixed-direction optimization in RTN components
Fixed-magnitude optimization (should be considered carefully because it is not ensured to find a solution to the problem if the thrust magnitude is too low)
Alternative PoC formulations (Chan's or Alfriend and Akella's)

The project is written in MATLAB and C++, and the required version of MATLAB is 2021b or newer. The optimization problem is run in MATLAB, whereas C++ functions are used to propagate the equations of motion of the spacecraft and compute the polynomial expansion of the probability of collision and other relevant quantities. The user is welcome to include their own dynamic models, but should do so only after understanding how the C++ and MATLAB environments interact. Therefore, it is advisable not to modify the C++ functions and only build the CAM scenarios using the MATLAB script. 
The MATLAB functions that should be modified to define the scenario are:
main.m
initOpt.m
defineParams.m

In main.m:
it is possible to define if the conjunctions are multiple and if they are in the cislunar domain by changing the two booleans in lines 3 and 4. 
In line 7, the firing times are given: these are "negative" times, meaning that they define the time that passes between the maneuver and TCA. If the conjunction is in Earth orbit, fireTimes are in orbital period units, whereas in cislunar they are in days. If more than one impulse is wanted, fireTimes should be a row vector including the times of each impulse: so, if, for example, we want to firing windows, one half an orbit before TCA and one 2 orbits after TCA (for a SK constraint), fireTimes = [0.5,-2]. If the selected thrust model is low-thrust, fireTimes should have at least two elements to determine the start and the end of the finite burn. If more than two elements are given, the algorithm detects if the firing times correspond to the same thrusting window (as is the case in Section 6.C of Ref. [1]) or if they define multiple thrust arcs. To define multiple thrust arcs, the time step between discretization nodes inside a thrust arc should be constant, and a larger time step should be given between the end time of one thrust arc and the start time of the subsequent one.

In initOpt.m
One can define the harmonics in line 35: gravOrd = 2, J2 dynamics;  gravOrd = 4, J4 dynamics. Also, one can insert the user-defined scenarios parameters of the conjunction: primary and secondary mean state and covariance at TCA, HBR, and physical properties if needed by the dynamics models. To do this one should create a "generateInit" function, using "generateInitLeo","generateInitMultiple" and "generateInitiCislunar" as examples.
The output of the "generateInit" function must be a structure with fields "mu","Lsc","Vsc","Tsc","n_conj","tca_sep","T","primary", and "secondary". mu is the gravitational parameter, Lsc, Vsc, and Tsc are scaling constants, n_conj is the number of consecutive conjunctions, tca_sep is the time difference between the first conjunction and all the subsequent ones, in seconds (first element must be 0), T is the orbital period in seconds (1 for Cislunar); "primary" and "secondary" are substructures wich include the state at TCA in ECI (synodic for Cislunar) x0, the position covariance matrix at TCA in ECI (synodic in Cislunar), the HBR of the body, and the physical parameters used for the envornmental models (if needed).

In defineParams.m
The parameters of the optimization are defined in lines 16-29: the order of the polynomial expansion, the type of PoC model approximation, the method to solve the NLP (recursive or fmincon), and the required constraint thresholds to be respected. Moreover, one can select the required operational constraints in lines 30-43.

To run the codes, the first thing is to download this repository and install the required Differential Algebra tool, DACE, which is publicly available at https://github.com/dacelib/dace. Afterwards, download this repository and unzip it to the preferred location. Compile the executables from validatePoly.cpp and polyProp.cpp via CMAKE. The CMakeLists.txt can be used to this end. The executables should be placed in the folder ./CppExec, so, after installing cmake in your bash, it is advisable to run "cmake .." and "make" inside this folder.

DOI: 10.5281/zenodo.11169701

Author: Zeno Pavanello, 2025

E-mail: zpav176@aucklanduni.ac.nz

[1] Pavanello, Z., Pirovano, L., & Armellin, R. (2025). Recursive Polynomial Method for Fast Collision Avoidance Maneuver Design. IEEE Transactions on Aerospace and Electronic Systems, 61(2), 2914–2925. https://doi.org/10.1109/TAES.2024.3480041
[2] Pavanello, Z., Pirovano, L., & Armellin, R. (2025). Efficient collision avoidance manoeuvres under multiple polynomial constraints. Acta Astronautica, 234, 548–559. https://doi.org/10.1016/j.actaastro.2025.05.002
