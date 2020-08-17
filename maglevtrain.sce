function y1=maglevtrain(u1,u2,u3,u4)
// This is nonlinear maglev trainmodel

// Load the parameters
// Maglev train parameters
A=0.024;
M=700;
N=450;
R=1.2
mu=4*%pi*(10^(-7));
g=9.81;
// Variables de estado
//u1;separacion
//u2;balanceo
//u3;densidad de flujo magnetico

// Variables de control
//u4;voltaje
//Ecuciones de variables de estado no lineal
y1=[u2;
       g-(((A)/(mu*M))*u3^2);
       -(((2*R)/(mu*A*N^2))*u1*u3)+ (u4/(N*A))]

endfunction 

