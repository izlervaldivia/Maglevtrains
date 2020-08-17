//Lqrtrain
clear
clc
A=0.024;
M=3;
N=450;
R=1.2
mu=4*%pi*(10^(-7));
g=9.81;
B=0.6
z=0.008
A=[0 1 0;0 0 -((2*B*A)/(mu*M));-((2*R*B)/(mu*A*(N^2))) 0 -((2*R*z)/(mu*A*(N^2)))]
B=[0;0;0.0925926];//[0;0;1/N*A]
C=[1 0 0
   0 0 1];
D=[0;0];
sys = syslin("c",A, B,C);
Q_xx=diag([100 1 1]); //Weights on states
R_uu   = 0.01; //Weight on input
G1=lqr(sys,Q_xx,R_uu);
