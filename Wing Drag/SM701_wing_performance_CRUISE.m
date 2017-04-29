% Aircraft drag build up computations
% Jack Langelaan, 2013

%{
You need two other functions to make this work:
ParseXFOILData_alpha
ComputeSingleSurfaceDrag
%}

clear all
close all
clc

d2r=pi/180;
r2d=180/pi;

% reference quantities
atmosphere.T=5+460; % R
atmosphere.rho=0.001496; % slug/ft^3
atmosphere.mu=0.000000343;% lb s/ft 
world.g=32.174; % ft/s^2;

% aircraft reference data
aircraft.S=120; % ft^2
aircraft.b=40; % ft
aircraft.m=70.24294424; % slug

% some reference quantities that are computed
aircraft.AR=aircraft.b^2/aircraft.S; % aspect ratio
aircraft.c=aircraft.S/aircraft.b; % mean geometric chord

% path to airfoil files
airfoil_name='SM701.dat';
wing_airfoil_file={
    'SM701/XFOIL/SM701_1e6.txt';
    'SM701/XFOIL/SM701_2e6.txt';
    'SM701/XFOIL/SM701_3e6.txt';
    'SM701/XFOIL/SM701_4e6.txt';
    'SM701/XFOIL/SM701_5e6.txt';
    'SM701/XFOIL/SM701_6e6.txt';
    'SM701/XFOIL/SM701_6e6.txt';%Use for Re# > 6e6
    };

% get airfoil data
alpha_range=-10:0.25:15; % range of alphas for which you have XFOIL data
CL_range=0:0.1:2; % range of CLs that you need to cover in your aircraft
wing_airfoil=ParseXFOILData_alpha(1,wing_airfoil_file,...
    [1e6 2e6 3e6 4e6 5e6 6e6 12e6],...
    alpha_range,CL_range,airfoil_name);


%CL cruise = .73
%CL max = 1.6
CL=0.2:0.1:1.6; % CL range. You will change this so it covers about 0.2 to 1.6 i.e. CL=0.2:0.1:1.6
%{
This next bit defines paths to your AVL output files. Make sure they are in
the same order as your CL vector in the line above!!!
%}
surface_files={
    'SM701/AVL/SM701_010';
    'SM701/AVL/SM701_020';
    'SM701/AVL/SM701_030';
    'SM701/AVL/SM701_040';
    'SM701/AVL/SM701_050';
    'SM701/AVL/SM701_060';
    'SM701/AVL/SM701_070';
    'SM701/AVL/SM701_080';
    'SM701/AVL/SM701_090';
    'SM701/AVL/SM701_100';
    'SM701/AVL/SM701_110';
    'SM701/AVL/SM701_120';
    'SM701/AVL/SM701_130';
    'SM701/AVL/SM701_140';
    'SM701/AVL/SM701_150';
    'SM701/AVL/SM701_160';
    'SM701/AVL/SM701_170';
    'SM701/AVL/SM701_180';
    };

NCL=length(CL);

DO_PLOTS=0; % this makes teh drag function plot lift coefficients. Set to 0 to suppress
for i=1:NCL
    surffile=surface_files{i};
    
    q(i)=aircraft.m*world.g/(aircraft.S*CL(i));
    v(i)=sqrt(2*q(i)/atmosphere.rho);
    vplot(i)=v(i)*0.592484;
    
    % average reynolds number
    Re(i)=atmosphere.rho*v(i)*aircraft.c/atmosphere.mu;
    
    % aerosurface drag
    [S_w,CD.prof(i),CD.ind(i)]=ComputeSingleSurfaceDrag(DO_PLOTS,atmosphere,CL(i),v(i),surffile,wing_airfoil);

    % You should include your fuselage drag computation here.
    D.fuse(i)=0; % MAKE SURE YOU REPLACE THIS WITH YOUR COMPUTED DRAG
    
    % Do your trim drag computation here.
    D.trim(i)=0; % MAKE SURE YOU REPLACE THIS WITH YOUR COMPUTED DRAG
    
    D.prof(i)=q(i)*aircraft.S*CD.prof(i);
    D.ind(i)=q(i)*aircraft.S*CD.ind(i);
    e(i)=(CL(i)^2/(pi*aircraft.AR))/CD.ind(i); % compute span efficiency
    
    D.tot(i)=D.ind(i)+D.prof(i)+D.fuse(i)+D.trim(i);

    CD.tot(i)=D.tot(i)/(q(i)*aircraft.S);    
end

fprintf('AVL computed wing area: %g\n',S_w)
    
figure('Name','Drag plot')
plot(vplot,D.ind,'r',vplot,D.prof,'g',vplot,D.fuse,'c',vplot,D.trim,'m',vplot,D.tot,'b')
legend('induced','profile','fuselage','trim','total');
xlabel('va (knots)')
ylabel('Drag force (lb)')

figure('Name','Drag Coefficient')
subplot(2,1,1)
plot(CL,CD.tot)
ylabel('C_D')
subplot(2,1,2)
plot(CL,CL./CD.tot)
ylabel('L/D')
xlabel('C_L')

 A=[ones(NCL,1) transpose(CL(1:NCL)).^2];
 CDpoly=(A'*A)\A'*CD.tot(1:(NCL))';
 fprintf('Second order fit: CD0= %g; Oswald efficiency= %g\n',CDpoly(1),1/(pi*aircraft.AR*CDpoly(2)))
 
