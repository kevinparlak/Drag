%{
Aircraft drag computation for AERSP402A
JWL 2015
%}
clc
close all
clear all

global world

SaveData=0; % set to 1 to save a lot of data to a file...
SaveSuffix='_afterdrop.mat'; % used for generating the save file name...

% some world parameters
%Takeoff parameters
world.g=32.174;         %ft/s^2
world.rho=0.001496;      %slug/ft^3
world.mu=0.000000343;   %lb s/ft^2

aircraft.name='SPYRO'; % used in save file name
%English mass empty
ENGm_empty=24.39854479;                     %slugs
aircraft.m_empty=ENGm_empty;                %slugs (includes structure, servos, wiring, ESC, autopilot)
%English mass payload
ENGm_payload=0;
aircraft.m_payload=ENGm_payload;            %slugs
%English mass battery
ENGm_battery=20.97967303;
aircraft.m_batt=ENGm_battery;               %slugs (battery)
aircraft.e_batt=260*3600*10.7639;           %lb ft/slug battery energy density
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
aircraft.S=ENGarea;                 %reference area, ft^2
aircraft.eta=0.8*0.9*0.9;           %drivetrain total efficiency (prop*motor*speedcontroller)
%Flat plate area
aircraft.fpar=.181994+.127595;           %0.013 is a square plate with area 0.01m (10cm x 10cm) x 3.28 ft/m perp. to flow
%English fuselage length
ENGlength=18;                       %ft
aircraft.fuse=ENGlength;            %ft fuselage length
%English fuselage wetted area
ENGwet=61.6+2*(27.5);                   %ft^2
aircraft.fuse_wetted=ENGwet;            %ft^2 fuse wetted area

savepath=strcat(aircraft.name,SaveSuffix);

% path to needed files
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
airfoil(1)=ParseXFOILData_alpha(1,front_airfoil_file,[3e5 1e6 2e6 3e6 4e6 5e6 6e6 15e6],-10:0.25:15,0:0.1:1.6,front_airfoil_name);
airfoil(2)=ParseXFOILData_alpha(1,rear_airfoil_file,[3e5 1e6 2e6 3e6 4e6 5e6 6e6 15e6],-10:0.25:15,-1:0.1:1.2,rear_airfoil_name);
airfoil(3)=ParseXFOILData_alpha(1,fin_airfoil_file,[3e5 1e6 2e6 3e6 4e6 5e6 6e6 15e6],-10:0.25:15,-1:0.1:1.2,fin_airfoil_name);


% path to strip force files. Notice that aircraft.name shows up here...
% that means that there is a directory called ../AVL/MyAircraft that
% contains all the AVL data, with strip force files living in
% MyAircraft/FS. I also assume that that the AVL files have all been saved 
% as ??.out, where ?? is 10*CL (so CL=0.4 is 04.out)

CL=0.3:0.1:1.6; % CL range
NCL=length(CL);

% some computed parameters
aircraft.c=aircraft.S/aircraft.b;
aircraft.AR=aircraft.b^2/aircraft.S;
aircraft.m_tot=aircraft.m_empty+aircraft.m_payload+aircraft.m_batt;         %total mass (slug)

forcefile={
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
    [f_ind(:,i),f_pro(:,i),f_lift(:,i),S_surf]=ComputeAeroSurfaceDrag_V2(DO_PLOTS,SurfaceDef,CL(i),v(i),surffile,airfoil);    
    D_ind(:,i)=q(i)*f_ind(:,i);
    D_pro(:,i)=q(i)*f_pro(:,i);
    L_surf(:,i)=q(i)*f_lift(:,i);
    
    fprintf('Computed lift force: %g; Load factor: %g. \n',sum(L_surf(1:4,i)),sum(L_surf(1:4,i))/(aircraft.m_tot*world.g));
    
    Re.fus(i)=world.rho*v(i)*(aircraft.fuse/2)/world.mu;
    cf=0.455*(log(Re.fus(i))^(-2.58));
    f_fus(i)=aircraft.fuse_wetted*cf;
    D_fus(i)=f_fus(i)*q(i);
    
    D_par(i)=q(i)*aircraft.fpar; % extra parasite drag
    
    D_tot(i)=sum(D_ind(:,i))+sum(D_pro(:,i))+D_fus(i)+D_par(i);
    
    aircraft.CL(i)=CL(i);
    aircraft.v(i)=v(i);
    aircraft.q(i)=q(i);
    aircraft.D_tot(i)=D_tot(i);
        
    P.aero(i)=D_tot(i)*v(i);
    P.batt(i)=P.aero(i)/aircraft.eta;
        
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

if SaveData
    save(savepath,'aircraft','world');
end


figure('Name','Drag Polar')
plot(CL,CD_tot,'ro',CL,polyval(CDpoly2,CL),'b')

% Now do some plotting.
v_u=v./1.6781;  %convert to knots
speedlabel='airspeed (knots)';
figure
plot(v_u,D_ind(1,:)+D_ind(2,:),'r',...
    v_u,D_ind(3,:)+D_ind(4,:),'g',...
    v_u,D_ind(5,:),'b',...
    v_u,D_pro(1,:)+D_pro(2,:),'r:',...
    v_u,D_pro(3,:)+D_pro(4,:),'g:',...
    v_u,D_pro(5,:),'b:',...
    v_u,D_fus,'c',...
    v_u,D_par,'m',...
    v_u,D_tot,'k');
legend('front wing induced','rear wing induced','tail induced',...
    'front wing profile','rear wing profile','tail profile',...
    'fuselage','parasite','total','Location','NorthWest')
ylabel('drag (lb)')
figure
plot(v_u,LD,'b')
ylabel('L/D');
xlabel(speedlabel)
figure
ha2=subplot(2,1,1);
plot(v_u,P.aero/550,'b')
ylabel('Aero power required (hp)')
subplot(2,1,2)
plot(v_u,LD,'b')
ylabel('L/D');
xlabel(speedlabel)


figure('Name','Drag and Power')
subplot(3,1,1)
plot(v,LD)
ylabel('L/D')
subplot(3,1,2)
plot(v,D_tot)
ylabel('Drag (lb)')
subplot(3,1,3)
plot(v,P.batt/550)
ylabel('power (hp)')
xlabel('airspeed (ft/s)')


