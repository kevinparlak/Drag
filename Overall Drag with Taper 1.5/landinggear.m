function D_landinggear=landinggear(V,rho,airfoil)
%-----Wheel----
%Cylinder Cd
Cd_cylinder=.3;
%Coefficient of drag from tire
Cd_tire=1.1;
%Fairing for wheel drag coefficient
Cd_tirefairing=.00725;
%Wheel fairing thickness (%)
wheel=24;

%-----Nose gear-----
%NACA0024
%Airfoilthickness (%)
airfoilthickness=12;
%Strut thickness
t_nose=4/12;             %ft
%Strut length
l_nose=34/12;            %ft
%Tire thickness
w_tire=6.1/12;          %ft
%Tire length
d_tire=14/12;           %ft
%Coefficient of drag from airfoil
Cd_nose=.006;    
%Flat plate area
if airfoil==1
    fnose=Cd_nose*l_nose*(t_nose/(airfoilthickness/100));
    ftire=Cd_tirefairing*d_tire*(w_tire/(wheel/100));
else
    fnose=Cd_cylinder*l_nose*t_nose;
    ftire=Cd_tire*d_tire*w_tire;
end

%-----Main gear-----
%Same airfoil
%Strut thickness
t_main=(3.5/12);      %ft
%Strut length
l_main=(31.5/12);     %ft
%Tire thickness
w_tiremain=(5.25/12);          %ft
%Tire length
d_tiremain=(12/12);           %ft
%Coefficient of drag from airfoil
Cd_main=.006;    
%Flat plate area
if airfoil==1
    fmain=2*Cd_main*l_main*(t_main/(airfoilthickness/100));
    ftiremain=2*Cd_tirefairing*d_tiremain*(w_tiremain/(wheel/100));
else
    fmain=2*Cd_cylinder*l_main*t_main;
    ftiremain=2*Cd_tire*d_tiremain*w_tiremain;
end

%Landing gear drag
D_landinggear=.5*rho*V^2*(fnose+ftire+fmain+ftiremain);
end
