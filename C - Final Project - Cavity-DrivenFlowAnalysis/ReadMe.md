Like all other projects, the Final Project was coded in C. First I had to develop the numerical schemes applicable to the problem presented to me. 
The procedure and the derivation are described in detail, in the word document named "Final Project Report". 

Once the equations were derived, I had to code and solve the cavity-driven inflow condition so the vorticity effect needed in the flow can be 
produced. This is also detailed in the project report document. It is a 2-D enclosure with delta x and delta y steps = 0.00625.
The velocity of the plate is given as 1 (units), with a Reynolds number of 1000. The Error Max < 0.001 is the given convergence criteria. 
The top plate is at 3350 pressure units. The pressure gradient at the left and the right plate is 0.
The simulation ran for 15,000 iterations with a time step of delta t = 0.003 (ND)
Line Gauss-Seidel was used for this simulation for the sake of simplicity, but at the price of computation time. 

Once the results were in, they were imported to MATLAB to create the streamline plots and the vortex plots presented. 
