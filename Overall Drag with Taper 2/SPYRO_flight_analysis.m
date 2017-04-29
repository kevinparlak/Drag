%----------------------------------SPYRO Flight Analysis----------------------------------
%Takeoff
%-----Parameters-----
g=32.174;               %ft/s^2
%-----SPYRO-----
%Max battery weight
W_batt_max=675;         %lb
%Max takeoff weight
w_tot=2260;             %lb
m_tot=2260/g;           %slug



P_runway=225;           %hp
V_runway=120;           %ft/s
D_runway=148;           %lb
%Propeller efficiency
etaprop=.85;
%Simulate runway takeoff
for i=0:6
E_i(i)=.5*m_tot*


if CL(i)>1
    fprintf('Runway cleared')
end
end

%Rate of climb
ROC=900;                %ft/s
D_climb
%Cruise
D_cruise
%Return
D_return
D_descend






%Plot of drag v time for the mission
figure
%Plot of airspeed v time for the mission
figure
%Plot of power v time for the mission
figure