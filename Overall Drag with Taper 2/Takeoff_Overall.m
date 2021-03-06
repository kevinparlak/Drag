%{
Aircraft drag computation for AERSP402A
JWL 2015
%}
clc
close all
clear all


global world

SaveData=0; % set to 1 to save a lot of data to a file...
SaveSuffix='_takeoff.mat'; % used for generating the save file name...

% some world parameters
%Takeoff parameters
world.g=32.174;         %ft/s^2
world.rho=.002377;      %slug/ft^3
world.mu=0.000000375;   %lb s/ft^2

aircraft.name='SPYRO'; % used in save file name

%English mass empty
ENGm_empty=24.39854479;                     %slugs
aircraft.m_empty=ENGm_empty;                %slugs (includes structure, servos, wiring, ESC, autopilot)
%English mass payload
ENGm_payload=24.86475902;
aircraft.m_payload=ENGm_payload;            %slugs
%English mass battery
ENGm_battery=20.97967303;
aircraft.m_batt=ENGm_battery;               %slugs (battery)
aircraft.w_batt=aircraft.m_batt*world.g;    %lb (battery)
%-----Battery energy density-----
energydensity=260;
aircraft.e_batt=energydensity*3600*10.7639;           %lb ft/slug battery energy density
aircraft.d_batt=0.75;                       %battery discharge fraction (full discharge is bad for the battery...)

SurfaceDef.N_surf=6; % 6 surfaces defined in AVL file (includes YDUP surfaces)
SurfaceDef.name={'RF_Wing','LF_Wing','RR_wing','LR_wing','VertTail','VertTail2'};
SurfaceDef.foil=[1 1 2 2 3 3]; % airfoil indices- defines which airfoil is used for each surface.

% aircraft reference parameters. These have to match the reference
% parameters in the AVL file.
%English wingspan
ENGwingspan=40;                     %ft
aircraft.b=ENGwingspan;             %reference wingspan, ft
%English area
ENGarea=120;                        %ft^2  
aircraft.S=ENGarea;                %reference area, ft^2
aircraft.eta=0.8*0.9*0.9;           %drivetrain total efficiency (prop*motor*speedcontroller)
%Flat plate area
aircraft.fpar=.181994+.127595;           %0.013 is a square plate with area 0.01m (10cm x 10cm) x 3.28 ft/m perp. to flow

%Takeoff distance
d=700;                              %ft
%Rolling resistance
cr=.05;
%Propeller efficiency
etaprop=.85;

%Airfoil on landing gear? ==1
foil=1;

savepath=strcat(aircraft.name,SaveSuffix);

%path to needed files
front_airfoil_name='SM701/XFOIL/SM701.dat';
front_airfoil_file={
    'SM701/XFOIL/SM701_1e6.txt';% use Re=1e6 polar for all Re < 1e6...
    'SM701/XFOIL/SM701_1e6.txt';
    'SM701/XFOIL/SM701_2e6.txt';
    'SM701/XFOIL/SM701_3e6.txt';
    'SM701/XFOIL/SM701_4e6.txt';
    'SM701/XFOIL/SM701_5e6.txt';
    'SM701/XFOIL/SM701_6e6.txt';
    'SM701/XFOIL/SM701_6e6.txt'; % use Re=6e6 polar for all Re > 6e6...
    };
rear_airfoil_name='NACA0009/XFOIL/NACA_0009.dat';
rear_airfoil_file={
    'NACA0009/XFOIL/NACA0009_1e6.txt';% use Re=1e6 polar for all Re < 1e6...
    'NACA0009/XFOIL/NACA0009_1e6.txt';
    'NACA0009/XFOIL/NACA0009_2e6.txt'; 
    'NACA0009/XFOIL/NACA0009_3e6.txt';
    'NACA0009/XFOIL/NACA0009_4e6.txt';
    'NACA0009/XFOIL/NACA0009_5e6.txt';
    'NACA0009/XFOIL/NACA0009_6e6.txt';
    'NACA0009/XFOIL/NACA0009_6e6.txt';% use Re=6e6 polar for all Re > 6e6...
    };
fin_airfoil_name='NACA0009/XFOIL/NACA_0009.dat';
fin_airfoil_file={
    'NACA0009/XFOIL/NACA0009_1e6.txt'; % use Re=1e6 polar for all Re < 1e6...
    'NACA0009/XFOIL/NACA0009_1e6.txt';
    'NACA0009/XFOIL/NACA0009_2e6.txt'; 
    'NACA0009/XFOIL/NACA0009_3e6.txt';
    'NACA0009/XFOIL/NACA0009_4e6.txt';
    'NACA0009/XFOIL/NACA0009_5e6.txt';
    'NACA0009/XFOIL/NACA0009_6e6.txt';
    'NACA0009/XFOIL/NACA0009_6e6.txt';% use Re=6e6 polar for all Re > 6e6...
    };

% get airfoil data
airfoil(1)=ParseXFOILData_alpha(1,front_airfoil_file,...
    [3e5 1e6 2e6 3e6 4e6 5e6 6e6 15e6],-10:0.25:15,0:0.1:1.6,front_airfoil_name);
airfoil(2)=ParseXFOILData_alpha(1,rear_airfoil_file,...
    [3e5 1e6 2e6 3e6 4e6 5e6 6e6 15e6],-10:0.25:15,-1:0.1:1.2,rear_airfoil_name);
airfoil(3)=ParseXFOILData_alpha(1,fin_airfoil_file,...
    [3e5 1e6 2e6 3e6 4e6 5e6 6e6 15e6],-10:0.25:15,-1:0.1:1.2,fin_airfoil_name);


% path to strip force files. Notice that aircraft.name shows up here...
% that means that there is a directory called ../AVL/MyAircraft that
% contains all the AVL data, with strip force files living in
% MyAircraft/FS. I also assume that that the AVL files have all been saved 
% as ??.out, where ?? is 10*CL (so CL=0.4 is 04.out)

CL=0.1:0.1:1.4; % CL range
NCL=length(CL);

% some computed parameters
aircraft.c=aircraft.S/aircraft.b;
aircraft.AR=aircraft.b^2/aircraft.S;
aircraft.m_tot=aircraft.m_empty+aircraft.m_payload+aircraft.m_batt;         %total mass (slug)
aircraft.w_tot=aircraft.m_tot*world.g;                                      %total weight (lb)

forcefile={
    'Wing and Tail AVL/fs_010';
    'Wing and Tail AVL/fs_020';
    'Wing and Tail AVL/fs_030';
    'Wing and Tail AVL/fs_040';
    'Wing and Tail AVL/fs_050';
    'Wing and Tail AVL/fs_060';
    'Wing and Tail AVL/fs_070';
    'Wing and Tail AVL/fs_080';
    'Wing and Tail AVL/fs_090';
    'Wing and Tail AVL/fs_100';
    'Wing and Tail AVL/fs_110';
    'Wing and Tail AVL/fs_120';
    'Wing and Tail AVL/fs_130';
    'Wing and Tail AVL/fs_140';
    'Wing and Tail AVL/fs_150';
    'Wing and Tail AVL/fs_160';
    'Wing and Tail AVL/fs_170';
    'Wing and Tail AVL/fs_180';
    };
    

DO_PLOTS=0; % set to 1 if you want to see plots of section lift coefficient
for i=1:NCL
    fprintf('*** Runnning CL=%g.\n',CL(i))
    surffile=strcat(forcefile{i});

    q(i)=aircraft.m_tot*world.g/(aircraft.S*CL(i));
    v(i)=sqrt(2*q(i)/world.rho);
    
    % average reynolds number
    Re.avg(i)=world.rho*v(i)*aircraft.c/world.mu;
    
    % aerosurface drag
    [f_ind(:,i),f_pro(:,i),f_lift(:,i),S_surf]=ComputeAeroSurfaceDrag_V2...
        (DO_PLOTS,SurfaceDef,CL(i),v(i),surffile,airfoil);    

    D_ind(:,i)=q(i)*abs(f_ind(:,i));
    D_pro(:,i)=q(i)*f_pro(:,i);
    L_surf(:,i)=q(i)*f_lift(:,i);
    
    fprintf('Computed lift force: %g; Load factor: %g. \n',sum(L_surf(1:4,i)),...
        sum(L_surf(1:4,i))/(aircraft.m_tot*world.g));
    
    %-----Other drag contributors-----
    %Fuselage and boom drag
    [f_fus(i),D_fus(i),f_boom(i),D_boom(i)]=skin_friction_drag(v(i),world.rho,world.mu);
    %Landing gear drag
    %Airfoil?  ==1
    D_land(i)=landinggear(v(i),world.rho,foil);
    %Extra parasite drag
    D_par(i)=q(i)*aircraft.fpar;
    
    %-----Total drag-----
    D_tot(i)=sum(D_ind(:,i))+sum(D_pro(:,i))+D_fus(i)+D_boom(i)+D_par(i)+D_land(i);
    aircraft.CL(i)=CL(i);
    aircraft.v(i)=v(i);
    aircraft.q(i)=q(i);
    aircraft.D_tot(i)=D_tot(i);

    
    %Power
    t(i)=d/v(i);                            %s
    %-----Energy-----
    %Rolling friction
    E_f(i)=cr*aircraft.m_tot*world.g;       %lb-ft
    %Drag
    E_d(i)=D_tot(i)*d;                      %lb-ft
    %Inertial energy
    E_i(i)=.5*aircraft.m_tot*v(i)^2;        %lb-ft
    %Total
    E(i)=(E_f(i)+E_d(i)+E_i(i))/etaprop;            %lb-ft
   
    %Required for takeoff
    P.req(i)=E(i)/t(i);                     %lb-ft/s
    
    %Aero power
    P.aero(i)=D_tot(i)*v(i);
    %Available
    P.batt(i)=(aircraft.m_batt*aircraft.e_batt*aircraft.d_batt)/t(i);
        
    LD(i)=aircraft.m_tot*world.g/D_tot(i);
    
    CD_tot(i)=D_tot(i)/(q(i)*aircraft.S);
    aircraft.CD_tot(i)=CD_tot(i);
    aircraft.LD(i)=LD(i);
    aircraft.P_aero(i)=P.aero(i);
    aircraft.D_tot(i)=D_tot(i);
end

% compute second order fit of drag data (to get Oswald's efficiency and
% CD0). Only look at the middle part of the data (ignore CL=1 and CL=0.9).
A=[transpose(CL(1:NCL-2)).^2 ones(NCL-2,1)];
CDpoly=(A'*A)\A'*CD_tot(1:(NCL-2))';
fprintf('Second order fit: CD0= %g; Oswald efficiency= %g\n',CDpoly(2),1/(pi*aircraft.AR*CDpoly(1)))

CDpoly2=polyfit(CL(1:NCL-2),CD_tot(1:(NCL-2)),2);

aircraft.CDpoly=CDpoly2;
aircraft.CD0=CDpoly(2); % parasite drag coefficient
aircraft.eO=1/(pi*aircraft.AR*CDpoly(1)); % Oswald's efficiency

% compute CL and airspeed for best L/D flight condition
aircraft.CL_nom=sqrt(aircraft.CD0*aircraft.AR*aircraft.eO*pi);
aircraft.v_nom=sqrt((2*aircraft.m_tot*world.g)/(world.rho*aircraft.S*aircraft.CL_nom));

%To get airborne in 700 ft, CL>1
error=1;
CLneeded=1.1;
Vneeded=sqrt((2*aircraft.w_tot)/(world.rho*aircraft.S*CLneeded));
index=find(abs(Vneeded-v)<=error);
Pneeded=P.req(index)/550;
Eneeded=E(index);
tneeded=t(index);
D=aircraft.D_tot(index);
LDbest=LD(index);

if SaveData
    save(savepath,'aircraft','world');
end


figure('Name','Drag Polar')
plot(CL,CD_tot,'ro',CL,polyval(CDpoly2,CL),'b')

% Now do some plotting.
v_u=v./1.6781;  %convert to knots
Vknotsneeded=v_u(index);
speedlabel='airspeed (knots)';
%Drag
figure
plot(v_u,D_ind(1,:)+D_ind(2,:),'r',...
    v_u,D_ind(3,:)+D_ind(4,:),'g',...
    v_u,D_ind(5,:),'b',...
    v_u,D_pro(1,:)+D_pro(2,:),'r:',...
    v_u,D_pro(3,:)+D_pro(4,:),'g:',...
    v_u,D_pro(5,:),'b:',...
    v_u,D_fus,'c',...
    v_u,D_boom,'k:',...
    v_u,D_land,'m:',...
    v_u,D_par,'m',...
    v_u,D_tot,'k',...
    Vknotsneeded,D,'*');
legend('front wing induced','rear wing induced','tail induced',...
    'front wing profile','rear wing profile','tail profile',...
    'fuselage','boom','landing gear','parasite','total','SPYRO takeoff drag',...
    'Location','NorthWest')
ylabel('drag (lb)')
xlabel('airspeed (knots)')
axis tight

%L/D v knots
figure
plot(v_u,LD,'b',Vknotsneeded,LDbest,'*')
ylabel('L/D');
xlabel(speedlabel)
legend('LD profile','SPYRO Takeoff')

%Aero power required (hp)
figure
ha2=subplot(2,1,1);
plot(v_u,P.aero/550,'b',Vknotsneeded,P.aero(index)/550,'*')
ylabel('Aero power required (hp)')
legend('Power req','SPYRO Takeoff')
subplot(2,1,2)
plot(v_u,LD,'b',Vknotsneeded,LDbest,'*')
ylabel('L/D');
xlabel(speedlabel)
legend('LD profile','SPYRO Takeoff')

%Power comparisons
figure('Name','Drag and Power')
subplot(3,1,1)
plot(v,LD,Vneeded,LDbest,'*')
ylabel('L/D')
legend('LD profile','SPYRO Takeoff')
subplot(3,1,2)
plot(v,D_tot,Vneeded,D,'*')
ylabel('Drag (lb)')
legend('Total drag','SPYRO Takeoff')
subplot(3,1,3)
plot(v,P.req/550,Vneeded,Pneeded,'*')
legend('Power req (hp)','SPRYO Takeoff')
ylabel('power (hp)')
xlabel('airspeed (ft/s)')
axis([100,150,0,300])

%Battery weight used for takeoff
W_batt_takeoff=(Eneeded/aircraft.e_batt)*world.g;       %lb

%-----Output-----
fprintf('\nPower to takeoff (hp) = %f\n',Pneeded)
fprintf('Airspeed to takeoff (ft/s) = %f\n',Vneeded)
fprintf('Time to takeoff (s) = %f\n',tneeded)
fprintf('Takeoff drag (lb) = %f\n',D)
fprintf('Battery used (lb) = %f\n',W_batt_takeoff)

%-----Analysis-----
timeplot=t(1:index);
%---Drag v time---
figure
suptitle('SPYRO Takeoff Analysis')
%Best drag
dragplot=D_tot(1:index);
subplot(2,2,1)
plot(timeplot,dragplot)
xlabel('Time (s)')
ylabel('Drag (lb)')
%---Battery weight v time---
%Battery weight used
batteryplot=aircraft.w_batt-(E(1:index)/aircraft.e_batt)*world.g; 
subplot(2,2,2)
plot(timeplot,batteryplot)
xlabel('Time (s)')
ylabel('Battery used (lb)')
%---Airspeed v time---
subplot(2,2,3)
%Velocity change in time
vplot=sqrt(2*(((P.req(1:index).*timeplot*etaprop))-E_f(1:index)-E_d(1:index))/aircraft.m_tot);
plot(timeplot,vplot)
xlabel('Time (s)')
ylabel('Airspeed (ft/s)')
%---Power v time---
subplot(2,2,4)
%Power required
Pplot=P.req(1:index);
plot(timeplot,Pplot/550)
xlabel('Time (s)')
ylabel('Power required (hp)')