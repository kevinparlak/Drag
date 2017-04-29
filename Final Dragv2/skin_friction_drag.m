function [f_fuse,D_fuse,f_boom,D_boom,Re_fuse]=skin_friction_drag(V,rho,mu)

%FUSELAGE
X_fuse=[-6,-5.316,-4.632,-3.947,-3.263,-2.579,-1.895,-1.211,...
    -0.526,0.158,0.842,1.526,2.211,2.895,3.579,4.263,4.947,5.632,6.316];

A_fuse=[0,1.173,3.273,4.635,4.806,4.806,4.798,4.788,4.785,4.785,...
    4.785,4.786,4.785,4.786,4.784,4.784,4.574,3.213,1.171];

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

A_boom=[0,0.753,1.672,2.29,2.352,2.221,2.078,1.983,1.856,1.73,1.637,...
    1.503,1.361,1.217,0.981,0.762,0.574,0.36,0.144];

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