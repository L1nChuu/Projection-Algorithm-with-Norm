# Introduction

## Summary
Based on the algorithms of Liu2020 and Thom2015, the code implementation of 1-norm and 2-norm ball or sphere intersection projection algorithms is given.  
*Liu2020:Projections onto the Intersection of a One-Norm Ball or Sphere and a Two-Norm Ball or Sphere*  
*Thom2015:Efficient Dictionary Learning with Sparseness-Enforcing Projections*
## Procedures
>1.Proj_B1B2.m  
This function is aimed to solve the optimization problem.
![P1Omega](https://github.com/L1nChuu/Projection-Algorithm-with-Norm/blob/master/Proj_B1B2.png)    

>2.Proj_S1S2.m  
This function is aimed to solve the optimization problem.
![P1Omega](https://github.com/L1nChuu/Projection-Algorithm-with-Norm/blob/master/Proj_S1S2.png)   

>3.Proj_B1S2.m
This function is aimed to solve the optimization problem.
![P1Omega](https://github.com/L1nChuu/Projection-Algorithm-with-Norm/blob/master/Proj_B1S2.png)   

>4.FindRoot_QASB.m   
Based on the Quadratic Approximation Secant Bisection Method in Liu2020, this function gives a way of solving the root of function Psi and projection problem.

>5.Psi_Thom.m   
This function corresponds to the algorithm 2 in Thom2015, which is a subfunction about solving projection problem with Thom's method. 

>6.FindRoot_Thom.m   
Corresponding to algorithm3 in Thom2015, this function is aimed to solve the root of function Psi and projection problem.
 
>7.eplb.mexw64  
