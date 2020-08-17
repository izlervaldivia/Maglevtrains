////Matrices obtenidas de Scilab///
// load the data
clc
clear
load("maglevtrainLTI.sod","X","U","sys")
A=sys.A
B=sys.B
C=sys.C
D=sys.D
/////Controllability and Observability/////
///Matrices Scilab///
//controlabilidad
Cc = cont_mat(A,B)
rankCc=rank(Cc)
//observabilidad
O = obsv_mat(A, C)
rankO=rank(O)
/////Plotear valores singulares//////
///Valores matrices Scilab///
G = syslin('c', A, B, C, D);
tr = trzeros(G)
w = logspace(-3,3);
sv = svplot(G,w);
//ploteo valores singulares//
scf(1);
plot2d("ln", w, 20*log(sv')/log(10))
xgrid(12)
xtitle("Singular values plot","Frequency (rad/s)", "Amplitude (dB)");
////Obtencion de la funcion de transferencia////
//MatricesScilab//
[h]=ss2tf(sys)
////Obtencion de polos y ceros/////
//Matrices Scilab//
scf(2)
plzr(h);
xtitle("Polos y zeros matrices Scilab")

