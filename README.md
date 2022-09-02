# LatticeOPT
A heuristics-based topology optimization algorithm for thin-walled lattice structures in Abaqus/Explicit

We include 3 examples for users to get familiarized with the use of LatticeOPT:
Example 1: Long, slender, lattice-filled column under axial compression
Example 2: Dynamic impact between a rigid pole and a lattice-filled beam fixed at both ends
Example 3: A lattice-filled sandwich panel subject to blast loading

In each example folder, we include a BoundaryCond.inp to set up the problem with appropriate boundary conditions.
The file InitialConditions.py sets up the initial condition and framework parameters for the example.

To run an example:
1. Move BoundaryCond.inp and InitialConditions.py to the same folder as the other source scripts
2. Run LatticeOPT in command line via: python LatticeOPT.py
3. Optimization results can be visualized by running PlotResults.py