function [f_fuse,D_fuse,f_boom,D_boom,Re_fuse]=skin_friction_drag(V,rho,mu)

%FUSELAGE
X_fuse=[-5,-4.026,-3.053,-2.079,-1.105,-0.132,0.842,1.816,2.789,3.763,...
    4.737,5.711,6.684,7.658,8.632,9.605,10.579,11.553,12.526];
A_fuse=[0,0.811,2.184,3.502,4.443,4.785,4.785,4.785,4.785,4.785,4.625,...
    4.28,3.962,3.506,2.926,2.424,1.653,0.93,0.31];
P_fuse=2.*sqrt(pi.*A_fuse);
Xle_fuse=X_fuse(1);
X_fuse=X_fuse-Xle_fuse;
a=length(X_fuse);

for j=1:a-1;
    dx_fuse(j)=X_fuse(j+1)-X_fuse(j);
    Swi_fuse(j)=(dx_fuse(j)./2)*(P_fuse(j)+P_fuse(j+1));
    Re_fuse(j)=(X_fuse(j))*V/mu;
    Cf_fuse(j)=0.455*(log10(Re_fuse(j)))^-2.58;
    f_fuse(j)=Swi_fuse(j)*Cf_fuse(j);
end
f_fuse=sum(f_fuse);
D_fuse=0.5*rho*V^2*f_fuse;

fprintf('f_fuse= %d (ft^2)\n',f_fuse)
fprintf('D_fuse= %d (lb)\n',D_fuse)

%BOOM
X_boom=[-1.5,-0.711,0.079,0.868,1.658,2.447,3.237,4.026,4.816,5.605,...
    6.395,7.184,7.974,8.763,9.553,10.342,11.132,11.921,12.711];

A_boom=[0,0.648,1.593,2.276,2.35,2.223,2.081,1.999,1.882,1.738,1.657,...
    1.532,1.372,1.248,1.008,0.763,0.622,0.398,0.158];
P_boom=2.*sqrt(pi.*A_boom);
Xle_boom=X_boom(1);
X_boom=X_boom-Xle_boom;
b=length(X_boom);

for j=1:b-1;
    dx_boom(j)=X_boom(j+1)-X_boom(j);
    Swi_boom(j)=(dx_boom(j)./2)*(P_boom(j)+P_boom(j+1));
    Re_boom(j)=(X_boom(j))*V/mu;
    Cf_boom(j)=0.455*(log10(Re_boom(j)))^-2.58;
    f_boom(j)=Swi_boom(j)*Cf_boom(j);
end
f_boom=sum(f_boom);
D_boom=0.5*rho*V^2*f_boom;

fprintf('f_boom= %d (ft^2)\n',f_boom)
fprintf('D_boom= %d (lb)\n',D_boom)
end