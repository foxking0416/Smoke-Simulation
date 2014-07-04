Daniel Rerucha and I are working as a team for our authoring tool project of CIS 660.
We use this smoke simulation framework as our based code for our authoring tool project.
We asked Dr. Kavan to get this smoke simulation framework a month before this assignment was assigned.
Therefore, our code would be similar.

Extra Credit:

1. Add preconditioner
	I implement the "Modified Incomplete Cholesky" preconditioner which is described in the fluids_notes.
	The function I modified is MACGrid::conjugateGradient.
	Without the preconditioner, the framerate is about 0.3 when simulating dimensions are 20 * 40 * 10
	With the preconditioner, the framerate is about 0.4 when simulating dimensions are 20 * 40 * 10	
2. Add cubic interpolation
	I implement the cubic interpolation method which is described in "Visual Simulation of Smoke" and fluids_notes.
	The function I modified is GridData::interpolate
3. Add fluid viscosity.
	I implement the fluid viscosity.
	I create a MACGrid::computeViscosityForce to compute the viscosity force and add the viscosity force in MACGrid::addExternalForces


Question
1.Describe the difference between Lagrangian and Eulerian viewpoints for simulation.  
  Why is the approach used in this assignment called Semi-Lagragian?
A: Lagrangian viewpoint:
       Quantities are carried by free moving particles. 
   Eulerian viewpoints:
       Quantities are meassured as it flows past.

	We called the method Semi-Lagragian because in order to compute the advect quantities, 
	we traced back the quantities of point which are originally located on the surface or center of each grid.
	The idea of tracing the quantities of a point is similar to the Lagrangian method, that's why we call it Semi-Lagrgian method.
  
2.Smoke and water can both be simulated with fluids.  
  Briefly explain how they are similar and how they are different. 
A: Smoke and water are defined as fluid, both flow following the pressure difference.
   They are different on the viscosity term and buoyancy force. 
   When we simulate the smoke, we usually discard the viscosity term and gravity term. 
   We also need to apply the buoyancy force on the smoke simulation.
   Instead, when we simulate water, we have to add the viscosity term and gravity term, and discard the buoyancy force.
   Besides, we also need to compute the surface tension of water surface.
  
3.List one advantage and one disadvantage to simulating our fluid on a grid.  
  Describe two other techniques for simulating fluids and the advantages and disadvantages of each.
A: 
   Eulerian grid method advantage: The ability of simulate the incompressibility property of fluid.
   Eulerian grid method disadvantage: Cannot simulate splashes, droplet.
   
   Lagrangian method advantage: Ability to simulate splashes and droplets.
   Lagrangian method disadvantage: Cannot simulate smooth surfaces.
   
   Lattice Boltzmann Methods
   
  
4.How do level set surfaces differ from mesh surfaces? 
  List one advantage and disadvantage of using a mesh to represent a surface.  
  List one advantage and disadvantage to using a Level Set Surface.
A: Mesh surface advantage: 
       Easy to compute and implement.
   Mesh surface disadvantage: 
       If we hope to have a smooth water surface, then the meshes are required to be super fine. 
	   That means crazy simulation time is required.
	   
   Level Set Surface advantage: 
       Very easy to model smooth surfaces. 
	   Provide information such as whether a point is inside water or not.
   Level Set Surface disadvantage:
       Cannot reliably represent small features. 
	   Level set represention tends to automatically eliminate small features as they move through the grid.