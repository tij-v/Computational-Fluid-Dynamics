## **Wave Behavior Simulation Inside a Closed Tube**

Task 1
Simulation in C to analyze shockwave propagation behavior inside a closed tube, using different numerical techniques. Namely:
1. Forward in Time Backward in Space (FTBS) Explicit Method
2. Lax-Wendroff Method
3. Euler's Backward in Time Center in Space (BTCS) Implicit Method
In this task, the simulation solely depends on the initial conditions and the boundary conditions. No particular governing model 
PDEs. As a subtask, I compared the effectiveness of each method at different time steps. The judging criteria were based on:
1. The loss of wave amplitude
2. The deformation of the wave shape over time.
The simulation ran for up to T = 0.15 seconds.

Task 2
Simulation in C to use MacCormack's Explicit Method to solve the Burger's Equation at different time steps:
1. delta t = 0.1
2. delta t = 0.2
The simulation ran for up to T = 2.4 seconds.
The results were plotted and then compared. 
 
