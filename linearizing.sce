// linearizing.sce
// Search the SUPERBLOCK in Xcos
for i=1:length(scs_m.objs)
    if typeof(scs_m.objs(i))=="Block" & scs_m.objs(i).gui=="SUPER_f" then
        scs_m = scs_m.objs(i).model.rpar;
        break;
    end
end

// Set the equilibrium point
X=[0.008;0;0.6];
U=[22.407];

// linearize the model
sys = lincos(scs_m,X,U);

// obtaingin the matrices A,B,C,D
A=sys.A 
B=sys.B
C=sys.C
D=sys.D

// save the data
save("maglevtrainLTI.sod","X","U","sys")
load("maglevtrainLTI.sod","X","U","sys")
