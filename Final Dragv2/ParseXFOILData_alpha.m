function airfoil=ParseXFOILData_alpha(SHOWPLOTS,filename,Re,alpha,cl_range,airfoil_name)
%{
Parse XFOIL data files to generate a matrix of cd vs cl and cm vs cl
This first interpolates over alpha because xfoil may have generated data
at extremely high alpha density (so that cl is constant to its precision)
and then interpolates over cl to give cd vs cl and cm vs cl.
Inputs:
SHOWPLOTS: 1= show plots, 0= don't plot
filename: set of full filenames (including paths) to the xfoil polar files
Re: vector of Reynolds numbers of the Polar files
alpha: the range of angles of attack over which you want to interpolate
data
cl_range: the range of cl over which you want airfoil data
airfoil_name: filename and path to the airfoi coordinates file (so the
airfoil can get plotted...

The output is a structure called airfoil:
airfoil.cl:     vector of cl
airfoil.Re:     vector of Re
airfoil.cd:     array of length(cl) x length(Re)
airfoil.cm:     array of length(cl) x length(Re)
airfoil.x       airfoil coordinates
airfoil.y
%}

airfoil.cl=cl_range;

airfoil.Re=Re;

for n=1:length(filename)
    fprintf('%s\n',filename{n}); % print the filename to the command window so you can see where you are
    j=0;
    fid=fopen(filename{n});
    for i=1:12, tline=fgets(fid); end % skip the first 12 lines
    while 1
        tline=fgets(fid);
        if tline==-1, break; end
        j=j+1;
        D=sscanf(tline,'%f');
        a(j)=D(1);
        cl(j)=D(2);
        cd(j)=D(3);
        cm(j)=D(5);
    end
    fclose(fid);
    a(j+1:length(a))=[];
    cl(j+1:length(cl))=[];
    cd(j+1:length(cd))=[];
    cm(j+1:length(cm))=[];
    cl_interp=interp1(a,cl,alpha,'spline',NaN);
    cd_interp=interp1(a,cd,alpha,'spline',NaN);
    cm_interp=interp1(a,cm,alpha,'spline',NaN);
    airfoil.cd(:,n)=interp1(cl_interp,cd_interp,airfoil.cl,'spline',NaN);
    airfoil.cm(:,n)=interp1(cl_interp,cm_interp,airfoil.cl,'spline',NaN);
    leg_text{n}=sprintf('Re= %0.2e',Re(n));
end

fid=fopen(airfoil_name);
airfoil.name=fgetl(fid);
j=0;
while 1
    tline=fgets(fid);
    if tline==-1, break; end
    j=j+1;
    D=sscanf(tline,'%f');
    airfoil.x(j)=D(1);
    airfoil.y(j)=D(2);    
end
fclose(fid);

if SHOWPLOTS
    figure('Name','Airfoil Drag Polar')
    subplot(2,1,1)
    plot(airfoil.x,airfoil.y)
    title(airfoil.name)
    axis equal
    subplot(2,1,2)
    plot(airfoil.cd,airfoil.cl)
    xlabel('c_d')
    ylabel('c_l')
    legend(leg_text);
end

return;