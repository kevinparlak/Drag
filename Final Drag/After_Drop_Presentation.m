%After-Drop Drag
clc
close all
clear all

global world

SaveData=0; % set to 1 to save a lot of data to a file...
SaveSuffix='_afterdrop.mat'; % used for generating the save file name...

%Atmospheric Parameters
world.g=32.174;         %ft/s^2
world.rho=0.001496;     %slug/ft^3
world.mu=0.000000343;   %lb-s/ft^2

aircraft.name='SPYRO'; % used in save file name

%-----MASS-----
%English mass empty
ENGm_empty=24.39854479;                     %slugs
aircraft.m_empty=ENGm_empty;                %slugs (includes structure, servos, wiring, ESC, autopilot)
%English mass payload
ENGm_payload=0;
aircraft.m_payload=ENGm_payload;            %slugs
%English mass battery
ENGm_battery=20.97967303;
aircraft.m_batt=ENGm_battery;               %slugs (battery)

%-----Battery energy density-----
energydensity=260;
aircraft.e_batt=energydensity*3600*10.7639; %lb-ft/slug battery energy density
aircraft.d_batt=0.75;                       %battery discharge fraction (full discharge is bad for the battery...)
energyvolume=350*75187.6;                   %battery volume (lb-ft/ft^3)

SurfaceDef.N_surf=6; % 6 surfaces defined in AVL file (includes YDUP surfaces)
SurfaceDef.name={'RF_Wing','LF_Wing','RR_wing','LR_wing','VertTail','VertTail2'};
SurfaceDef.foil=[1 1 2 2 3 3]; % airfoil indices- defines which airfoil is used for each surface.

%-----PARAMETERS-----
%English wingspan
ENGwingspan=40;                     %ft
aircraft.b=ENGwingspan;             %reference wingspan, ft
%English area
ENGarea=120.25;                     %ft^2
aircraft.S=ENGarea;                 %reference area, ft^2
aircraft.eta=0.85*0.9*0.9;          %drivetrain total efficiency (prop*motor*speedcontroller)
%Flat plate area
aircraft.fpar=.181994+.127595;      %0.013 is a square plate with area 0.01m (10cm x 10cm) x 3.28 ft/m perp. to flow

savepath=strcat(aircraft.name,SaveSuffix);

%Range
rangemiles=125;                     %miles
range=rangemiles*5280;              %ft

%Propeller efficiency
etaprop=.85;

%Retractable landing gear? ==1
gear=0;

%Airfoil on landing gear? ==1
foil=1;

%-----PATHS-----
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

%-----AIRFOIL DATA-----
airfoil(1)=ParseXFOILData_alpha(1,front_airfoil_file,...
    [3e5 1e6 2e6 3e6 4e6 5e6 6e6 15e6],-10:0.25:15,0:0.1:1.6,front_airfoil_name);
airfoil(2)=ParseXFOILData_alpha(1,rear_airfoil_file,...
    [3e5 1e6 2e6 3e6 4e6 5e6 6e6 15e6],-10:0.25:15,-1:0.1:1.2,rear_airfoil_name);
airfoil(3)=ParseXFOILData_alpha(1,fin_airfoil_file,...
    [3e5 1e6 2e6 3e6 4e6 5e6 6e6 15e6],-10:0.25:15,-1:0.1:1.2,fin_airfoil_name);

%-----STRIP FORCE FILES-----
%CL Range
CL=0.1:0.1:1.4;
NCL=length(CL);

%Computed Parameters
aircraft.c=1.13;
aircraft.AR=aircraft.b^2/aircraft.S;
aircraft.m_tot=aircraft.m_empty+aircraft.m_payload+aircraft.m_batt;         %total mass (slug)
aircraft.w_tot=aircraft.m_tot*world.g;                                      %total weight (lb)

%-----FORCE FILES-----
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
    
%-----PLOTS-----
DO_PLOTS=0; % set to 1 if you want to see plots of section lift coefficient
for i=1:NCL
    fprintf('*** Runnning CL=%g.\n',CL(i))
    surffile=strcat(forcefile{i});

    %Dynamic Pressure
    q(i)=aircraft.m_tot*world.g/(aircraft.S*CL(i));
    %Airspeed
    v(i)=sqrt(2*q(i)/world.rho);
    
    %Average Reynolds Number
    Re.avg(i)=world.rho*v(i)*aircraft.c/world.mu;
    
    %Aerosurface Drag
    [f_ind(:,i),f_pro(:,i),f_lift(:,i),S_surf]=ComputeAeroSurfaceDrag_V2...
        (DO_PLOTS,SurfaceDef,CL(i),v(i),surffile,airfoil);    
    %Induced Drag
    D_ind(:,i)=q(i)*abs(f_ind(:,i));
    %Profile Drag
    D_pro(:,i)=q(i)*f_pro(:,i);
    L_surf(:,i)=q(i)*f_lift(:,i);
    
    fprintf('Computed lift force: %g; Load factor: %g. \n',sum(L_surf(1:4,i)),...
        sum(L_surf(1:4,i))/(aircraft.m_tot*world.g));
    
    %-----DRAG CONTRIBUTORS-----
    %Fuselage and boom drag
    [f_fus(i),D_fus(i),f_boom(i),D_boom(i)]=skin_friction_drag(v(i),world.rho,world.mu);
    %Landing gear drag
    %Retractable? ==1
    if gear==1
        D_land(i)=0;
    else
        %Airfoil?  ==1
        D_land(i)=landinggear(v(i),world.rho,foil);
    end
    %Extra parasite drag
    D_par(i)=q(i)*aircraft.fpar;
    
    %-----TOTAL DRAG-----
    D_tot(i)=sum(D_ind(:,i))+sum(D_pro(:,i))+D_fus(i)+D_boom(i)+D_par(i)+D_land(i);
    %Aircraft Structure
    aircraft.CL(i)=CL(i);
    aircraft.v(i)=v(i);
    aircraft.q(i)=q(i);
    aircraft.D_tot(i)=D_tot(i);
    
    %-----POWER-----
    %Required Power
    P.aero(i)=D_tot(i)*v(i);
    P.aerohp(i)=P.aero(i)/550;
    %LD Ratio
    LD(i)=aircraft.m_tot*world.g/D_tot(i);
    CD_tot(i)=D_tot(i)/(q(i)*aircraft.S);
    
    aircraft.CD_tot(i)=CD_tot(i);
    aircraft.LD(i)=LD(i);
    aircraft.P_aero(i)=P.aero(i);
    aircraft.D_tot(i)=D_tot(i);
end

%-----MINIMUM DRAG-----
Dbest=min(D_tot);                               %lb
indexbest=find(D_tot==Dbest);
Vbest=v(indexbest);                             %ft/s
Pbest=P.aero(indexbest);                        %lb-ft/s
Ebest=Dbest*range/etaprop;                      %lb-ft
LDbest=LD(indexbest);
Wbest=(Ebest/aircraft.e_batt)*world.g;          %lb

%-----CONSTANT AIRSPEED FROM CRUISE MINIMUM DRAG-----
indexconstant=find(abs(v-187)<=10);
Dconstant=D_tot(indexconstant);                 %lb
Vconstant=v(indexconstant);                     %ft/s
Pconstant=P.aero(indexconstant);                %lb-ft/s
Econstant=Dconstant*range/etaprop;              %lb-ft
LDconstant=LD(indexconstant);
Wconstant=(Econstant/aircraft.e_batt)*world.g;  %lb

%Compute second order fit of drag data (to get Oswald's efficiency and
% CD0). Only look at the middle part of the data (ignore CL=1 and CL=0.9).
A=[transpose(CL(1:NCL-2)).^2 ones(NCL-2,1)];
CDpoly=(A'*A)\A'*CD_tot(1:(NCL-2))';
fprintf('Second order fit: CD0= %g; Oswald efficiency= %g\n',CDpoly(2),1/(pi*aircraft.AR*CDpoly(1)))

CDpoly2=polyfit(CL(1:NCL-2),CD_tot(1:(NCL-2)),2);

aircraft.CDpoly=CDpoly2;
%Parasite Drag Coefficient
aircraft.CD0=CDpoly(2);
%Oswald's Efficiency
aircraft.eO=1/(pi*aircraft.AR*CDpoly(1));

%CL and airspeed for best L/D flight condition
aircraft.CL_nom=sqrt(aircraft.CD0*aircraft.AR*aircraft.eO*pi);
aircraft.v_nom=sqrt((2*aircraft.m_tot*world.g)/(world.rho*aircraft.S*aircraft.CL_nom));

if SaveData
    save(savepath,'aircraft','world');
end

%-----DRAG POLAR-----
figure('Name','Drag Polar')
plot(CL,CD_tot,'ro',CL,polyval(CDpoly2,CL),'b')

%-----PLOT-----
v_u=v./1.6781;  %convert to knots
Vbestknots=v_u(indexbest);
Vconstantknots=v_u(indexconstant);
speedlabel='airspeed (knots)';

%-----DRAG-----
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
    Vbestknots,Dbest,'*',...
    Vconstantknots,Dconstant,'o');
legend('front wing induced','rear wing induced','tail induced',...
    'front wing profile','rear wing profile','tail profile',...
    'fuselage','boom','landing gear','parasite','total','After-Drop','Constant Airspeed from Cruise','Location','NorthWest')
ylabel('Drag (lb)')
xlabel('Airspeed (knots)')
axis tight

%-----LD-----
figure
plot(v_u,LD,'b',Vbestknots,LDbest,'*',Vconstantknots,LDconstant,'o')
ylabel('L/D');
xlabel(speedlabel)
legend('LD profile','Minimum Drag','Constant Airspeed from Cruise')

%-----AERODYNAMIC POWER REQUIRED-----
figure
ha2=subplot(2,1,1);
plot(v_u,P.aero/550,'b',Vbestknots,P.aero(indexbest)/550,'*',Vconstantknots,P.aero(indexconstant)/550','o')
ylabel('Aero power required (hp)')
legend('Power required (hp)','Minimum Drag','Constant Airspeed from Cruise')
subplot(2,1,2)
plot(v_u,LD,'b',Vbestknots,LDbest,'*',Vconstantknots,LDconstant,'o')
ylabel('L/D');
xlabel(speedlabel)
legend('LD profile','Minimum Drag','Constant Airspeed from Cruise')

%-----POWER COMPARISONS-----
figure('Name','Drag and Power')
subplot(3,1,1)
plot(v,LD,Vbest,LDbest,'*',Vconstant,LDconstant,'o')
ylabel('L/D')
legend('LD profile','Minimum Drag','Constant Airspeed from Cruise')
subplot(3,1,2)
plot(v,D_tot,Vbest,Dbest,'*',Vconstant,Dconstant,'o')
ylabel('Drag (lb)')
legend('Total drag','Minimum Drag','Constant Airspeed from Cruise')
subplot(3,1,3)
plot(v,P.aero/550,Vbest,Pbest/550,'*',Vconstant,Pconstant/550,'o')
legend('Power required (hp)','Minimum Drag','Constant Airspeed from Cruise')
ylabel('Power (hp)')
xlabel('Airspeed (ft/s)')
axis([100,400,0,500])

%-----Output-----
fprintf('\n\nAfter-Drop Analysis...')
%Minimum Drag
fprintf('\nMinimum Drag...')
fprintf('\nPower after-drop (hp) = %f\n',Pbest/550)
fprintf('Airspeed after-drop (ft/s) = %f\n',Vbest)
fprintf('After-drop drag (lb) = %f\n',Dbest)
fprintf('LD = %f\n',LDbest)
fprintf('Battery used (lb) = %f\n',Wbest)
fprintf('Range (miles) = %f\n',rangemiles)
%Constant Airspeed from Cruise
fprintf('\nConstant Airspeed from Cruise...')
fprintf('\nPower after-drop (hp) = %f\n',Pconstant/550)
fprintf('Airspeed after-drop (ft/s) = %f\n',Vconstant)
fprintf('After-drop drag (lb) = %f\n',Dconstant)
fprintf('LD = %f\n',LDconstant)
fprintf('Battery used (lb) = %f\n',Wconstant)
fprintf('Range (miles) = %f\n',rangemiles)