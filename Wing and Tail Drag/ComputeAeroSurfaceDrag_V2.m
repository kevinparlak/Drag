function [f_ind,f_pro,f_lift,S_surf]=ComputeAeroSurfaceDrag_V2(SHOWPLOTS,SurfaceDef,CL_nom,va,surfacefile,airfoil)
%{
Compute aerodynamic surface drag. This assumes symmetric flight (i.e. no sideslip)

coded:
April 23, 2015
Jack W. Langelaan

Inputs
SHOWPLOTS    show some plots

SurfaceDef
    .N_surf    number of surfaces (this includes duplicate surfaces!!)
    .name surface names (cell array)
    .foil airfoil index (array 1xN)

va            airspeed

surfacefile   name of strip force file

airfoil(N)    array of airfoil structures
    .Re    Reynolds number
    .cl    section cl
    .cd    section cd


Output (returns the drag area for each lifting surface)
f_ind    induced drag area array (1 x N)
f_pro    profile drag area array (1 x N)
f_lift   lift force area array (1 x N) (this is cl*s_surf)

%}

global world

nu=world.mu/world.rho; % kinematic viscosity

fid=fopen(surfacefile);

if SHOWPLOTS
    figure('Name',sprintf('Lift distributions: CL_nom=%d ',CL_nom))
    plothandle(1)=subplot(3,1,1);
    ylabel('c_l')
    hold on
    plothandle(2)=subplot(3,1,2);
    ylabel('c*c_l')
    hold on
    plothandle(3)=subplot(3,1,3);
    ylabel('Re (millions)')
    xlabel('y')
    hold on
end

f_pro=zeros(1,SurfaceDef.N_surf);
f_ind=zeros(1,SurfaceDef.N_surf);
f_lift=zeros(1,SurfaceDef.N_surf);
S_surf=zeros(1,SurfaceDef.N_surf);
for ns=1:SurfaceDef.N_surf
    foil_ID=SurfaceDef.foil(ns); % airfoil ID
    surfID=''; % find the next surface
    while (strcmp('  Surface #',surfID)==0)
        tline=fgets(fid);
        if length(tline)<14
            surfID='';
        else
            surfID=tline(1:11);
        end
    end
    fprintf('Analyzing %s',tline)
    fprintf('   Using airfoil: %s\n',airfoil(foil_ID).name)
    
    tline=fgets(fid);
    NS=sscanf(tline(37:44),'%d'); % number of sections in this surface
    
    for i=1:12, fgets(fid); end % skip 12 lines
    
    y=zeros(1,NS);
    c=zeros(1,NS);
    A=zeros(1,NS);
    cl=zeros(1,NS);
    cd_i=zeros(1,NS);
    Re=zeros(1,NS);
    cd_p=zeros(1,NS);
    for i=1:NS
        tline=fgets(fid);
        D=sscanf(tline,'%f');
        y(i)=D(2);
        c(i)=D(3);
        A(i)=D(4);
        cl(i)=D(8);
        cd_i(i)=D(9);
        Re(i)=va*c(i)/nu;
        cd_p(i)=interp2(airfoil(foil_ID).Re,airfoil(foil_ID).cl,airfoil(foil_ID).cd,Re(i),cl(i));        
    end
    S_surf(ns)=sum(A); % compute surface area
    f_pro(ns)=sum(A.*cd_p); % profile drag area
    f_ind(ns)=sum(A.*cd_i); % induced drag area
    f_lift(ns)=sum(A.*cl);  % lift area
    
    if SHOWPLOTS
        plot(plothandle(1),y,cl,'b')
        plot(plothandle(2),y,c.*cl,'b')
        plot(plothandle(3),y,Re/1e6,'b')
    end

end



fclose(fid);



return;
