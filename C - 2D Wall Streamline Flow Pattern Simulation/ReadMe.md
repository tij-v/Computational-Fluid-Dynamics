Simulation in C to capture the stream inflow behavior, give the wall conditions. 
This program implements the 2-D, incompressible flow, continuity partial differential equation to solve the problem. Discretization of this equation is done by implmenting various computational numerical schemes:
1. Line Gauss-Seidel Method
2. Point Gauss-Seidel Method
3. Line Successive Over Relaxation Method (LSOR)
4. Point Successive Over Relaxation Method (PSOR)

As a subtask, I also compared how soon each of the aforementioned methods solve the problem. This represented in the number of iterations given on each corresponding plot.
Furthermore, I also visualized at which iteration count does the program diverges and gets stuck in an infinite loop. Also given by a plot represenation. 
The simulation ran until the convergence criteria of Max Error < 0.01 was met.

