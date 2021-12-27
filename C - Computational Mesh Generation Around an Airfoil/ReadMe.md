## **Computational Elliptical Mesh Generation Around an Airfoil**

This is thus far my favorite coded simulation. For this simulation, not did we have to code an airfoil shape, but we also had to code the grid/mesh around it. 
This was particularly challenging because I had to connect all the points on the airfoil to all the points on the outer grid. 
This is plotted as the Algebraic Grid mesh. 
This was just the setup (initial condition).

Once the setup was ready, I was tasked to code and apply the Line Gauss-Seidel or Line Successive Over Relaxation Methods to transform the current Algebraic Mesh into an 
Elliptic Mesh type. The advantage of this mesh type is an increase in the accuracy of the calculations and a decrease in the computing time. 
This was plotted as the Elliptical Grid Mesh, where when compared to its counterpart, the difference is obvious. 
Line Gauss-Seidel was applied for the sake of simplicity and the simulation ran until the Max Error < 0.01 criteria was satisfied. 



