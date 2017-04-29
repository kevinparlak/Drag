function [f_fuse,D_fuse,f_boom,D_boom,Re_fuse]=skin_friction_drag(V,rho,mu)

%FUSELAGE
X_fuse=[-5,-4.0530,-3.1050,-2.1580,-1.2110,-0.2630,0.6840,1.6320,2.5790,3.5260,...
    4.4740,5.4210,6.3680,7.3160,8.2630,9.2110,10.1580,11.1050,12.0530];
A_fuse=[0,0.5220,1.1210,1.5010,1.6800,1.7300,1.7610,1.7630,1.6760,1.4080,1.1490,...
    1.1420,1.1500,1.1550,1.1250,0.7860,0.5230,0.3410,0.1640];
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
X_boom=[-1,-0.2630,0.4740,1.2100,1.9470,2.6840,3.4210,4.1580,4.8950,5.6320,...
    6.3680,7.1050,7.8420,8.5790,9.3160,10.0530,10.7900,11.5260,12.2630];

A_boom=[0,0.6840,1.2860,1.5460,1.5690,1.5110,1.4060,1.2260,1.0360,0.8990,0.8320,...
    0.7450,0.6200,0.4920,0.3980,0.3860,0.3750,0.3160,0.1760];
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