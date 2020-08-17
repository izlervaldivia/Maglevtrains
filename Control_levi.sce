//Evaluacion LQG Train
// load the data
clc
clear
load("maglevtrainLTI.sod","X","U","sys")

Ap=sys.A;
Bp=sys.B;
Cp=sys.C;
Dp=sys.D;
Dp=0
Cp=[1 0 0]
tri = trzeros(sys)
w = logspace(-3,3);
svi = svplot(sys,w);
scf(1);
plot2d("ln", w, 20*log(svi')/log(10))
xgrid(12)
xtitle("Valores singulares de la planta inicial","Frequency (rad/s)", "Amplitude (dB)");
///////////////////////////////////---------------------------

//Planta aumentada con el integrador
[ns,nc]=size(Bp); //ns= number of inputs; nc=number of controls
Ai=[Ap             Bp;
    0*ones(nc,ns) 0*ones(nc,nc)];

Bi=[0*ones(ns,nc); eye(nc)];
    
Ci=[Cp 0*ones(1,1)];

Di=0*ones(nc,nc);

sysi=syslin('c',Ai,Bi,Ci,Di);

I=eye(nc);

/*   Plot singular values    */
tri = trzeros(sysi)
w = logspace(-3,3);
svi = svplot(sysi,w);
scf(2);
plot2d("ln", w, 20*log(svi')/log(10))
xgrid(12)
xtitle("Valores singulares de la planta con integrador","Frequency (rad/s)", "Amplitude (dB)");
//Obtenciion de los polos y zeros planta con integrador
scf(3);
plzr(sysi);

//----LQR------//
//we use ricatti equation for calculate de gain H
C=1*Ci'*Ci;        //State Weighting Matrix
rho=1;       //Cheap control recovery parameter 
                //The smaller the parameter, the better the recovery.
R = rho*eye(nc);//Control Weigthing Matrix
//now we calculate B
B=Bi*inv(R)*Bi';
A=Ai;
//Solv the ricatti equation
X=riccati(A,B,C,'c','eigen');
//the value of the gain G
G=inv(R)*Bi'*X; //Matriz G
//----KALMAN FILTER-------///
ll= inv(Cp*inv(-Ap)*Bp+Dp);      //Choose ll and lh to match singular values at all frequencies
lh = -inv(Ap)*Bp*ll;

Lp=[lh,
   ll];                    //ll, lh - for low and high frequency loop shaping
   
pnint = eye(nc,nc)      // Process Noise Intensity Matrix
mu = 0.1;           // Measurement Noise Intesity; Used to adjust Kalman Filter Bandwidth
                     //Small mu - expensive sensor   - large bandwidth
                     //Large mu - inexpensive sensor - small bandwidth
THETA = mu*eye(nc,nc)   // Measurement Noise Intensity Matrix 

//We use the ricatti equation for calculate de gain H
Ch=Lp*Lp';
Ah=Ai';
//calculating Bh
Bh=Ci'*inv(THETA)*Ci;

//Calculate de solution
Xh=riccati(Ah,Bh,Ch,'c','eigen');

//The gain H
H=(inv(THETA)*Ci*Xh)';

sysh = syslin('c',Ai,H,Ci,Di);

/* Plot singular values*/
trh = trzeros(sysh)
w = logspace(-3,3);
svh = svplot(sysh,w);
scf(4);
plot2d("ln", w, 20*log(svh')/log(10))
xgrid(12)
xtitle("Valores singulares Malla objetivo G_{KF}", "Amplitude (dB)");
//--------------------------------------
//Compensator LQG
Ak = [ Ai-Bi*G-H*Ci  0*ones(ns+nc,nc)
       G          0*ones(nc,nc) ]

Bk = [ H
       0*ones(nc,nc) ]

Ck = [0*ones(nc, ns+nc) eye(nc,nc) ]

Dk = 0*ones(nc,nc);

sysk=syslin('c',Ak,Bk,Ck,Dk);

/* Plot singular values  */
trk = trzeros(sysk)
w = logspace(-3,3);
svk = svplot(sysk,w);
scf(5);
plot2d("ln", w, 20*log(svk')/log(10))
xgrid(12)
xtitle("Valores singulares compensador","Frequency (rad/s)", "Amplitude (dB)");


//----------------------------------------
//Analysis in open loop
Aol = [ Ap                     Bp*Ck
       0*ones(ns+nc+nc,ns)    Ak    ]

Bol = [ 0*ones(ns,nc)
       Bk ]
    
Col = [ Cp  0*ones(nc,ns+nc+nc) ]

Dol = 0*ones(nc,nc);

sysol = syslin('c',Aol,Bol,Col,Dol);

//Obtencion de los polos y zeros lazo abierto
scf(6);
plzr(sysol);
//----------------------------------------
//Response in closed loop
syscl = syslin('c',Aol-Bol*Col, Bol*0.0001, Col, 0*eye(nc,nc));
//Obtencion de los polos y  zeros en lazo cerrado
scf(7);
plzr(syscl);
//--------------------------------
//Respuesta al step
t=[0:0.1:20];
//input defined by a time function
deff('u=timefun(t)','u=1')
scf(8);
plot2d(t',(csim(timefun,t,syscl))')
xtitle("Respuesta del sistema","t(s)","Amplitud(m)");

