# Introduction

## Summary
Based on the Quadratic Approximation Secant Bisection  (QASB) Method in Liu2020 and the Bisection-Newton algorithm in Thom2015, the code implementation of the projection subproblem on the intersection of l1-norm ball/ sphere and l2-norm ball/sphere is given.
*Liu2020:Projections onto the Intersection of a One-Norm Ball or Sphere and a Two-Norm Ball or Sphere*  
*Thom2015:Efficient Dictionary Learning with Sparseness-Enforcing Projections*
## Procedures
>1.Proj_B1B2.m  
This function aims to solve the following projection suproblem for a given vector z. 
![P1Omega](https://github.com/L1nChuu/Projection-Algorithm-with-Norm/blob/master/Proj_S1S2.png)    

>2.Proj_S1S2.m  
This function aims to solve the following projection suproblem for a given vector z. 
![P1Omega](https://github.com/L1nChuu/Projection-Algorithm-with-Norm/blob/master/Proj_B1B2.png)   

>3.Proj_B1S2.m
This function aims to solve the following projection suproblem for a given vector z. 
![P1Omega](https://github.com/L1nChuu/Projection-Algorithm-with-Norm/blob/master/Proj_B1S2.png)   

>4.FindRoot_QASB.m   
Find a root of the function phi defined  in Liu20 using the Quadratic Approximation Secant Bisection Method (Algorithm 1) in Liu20.

>5.Psi_Thom.m   
This function corresponds to the algorithm 2 in Thom2015, which is a subfunction about solving projection problem with Thom's method. 

>6.FindRoot_Thom.m   
This function aims to find a  root of the function Psi defined in Thom2015, and solve the projection on the intersetion of l1-sphere and l2-sphere corresponding to algorithm3 in Thom2015.
 
>7.eplb.mexw64  
