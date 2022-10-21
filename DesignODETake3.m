function [Homework8_ODE] = DesignODETake3(V,Y)

%% Reactor Volume Parameters
% this comes first because everything is going to be done on a per tube basis
numTubes = 10; %CHANGE ME
tubeDiameter = 1; %meters, CHANGE ME

%% Other things that change
T = 533; %Kelvin
T_0 = Y(9); %Kelvin
P = 541.0755; %kPa
P_0 = 2000; %kPa

F_Prod = Y(1)/numTubes; %product
F_C2H4 = Y(2)/numTubes; %
F_HCl = Y(3)/numTubes; %
F_O2 = Y(4)/numTubes; %
F_CO2= Y(5)/numTubes;
F_H2O =Y(6)/numTubes;
F_Cl3Eth = Y(7)/numTubes;
F_Cl2 =Y(8)/numTubes;

F_tot = F_Prod + F_C2H4 + F_HCl + F_O2 + F_CO2 + F_H2O + F_Cl3Eth + F_Cl2; 

%% Ergun Related Equations

F_tot_0_overall = 970240 + 2910750 + 970240 + 5820; %[mol/hr], constant
F_tot_0 = F_tot_0_overall/numTubes; %[mol/hr] per tube, constant

rho_0 = 1.002; %kg/m3, constant
rho = rho_0*(P/P_0)*(T_0/T)*(F_tot_0/F_tot); %kg/m3, changing

volumetricFlowRate_tot_0 = 166430.8; %m3/hr, calculated in writeup, overall, constant
volumetricFlowRate_0 = volumetricFlowRate_tot_0/numTubes; %m3/hr, per tube, constant

Tau = 3.6; %space time in seconds, given 
reactorVol = Tau/3600*volumetricFlowRate_tot_0; %3600 to convert from s to hr
Ac=pi*(tubeDiameter/2)^2; %cross sectional area
%here, Matlab calculates a reactor length
% intentionally without a semicolon so it comes out in our command window
% as a result
reactorLength = reactorVol/numTubes/Ac;

superFicVelocity = volumetricFlowRate_0/Ac; % m/hr
G = rho*superFicVelocity; % superficial mass velocity [kg/m2-hr]

phi = 0.50; %void fraction/porosity [unitless], constant
mu = 0.023; %viscosity of gas [mPa-s], constant
gc = 1; %unit conversion [metric], constant
particleDiameter = tubeDiameter/8; %based on heuristic [m], constant

beta_0 = G*(1-phi)/(rho_0*gc*particleDiameter*phi^3)*(150*(1-phi)*mu/particleDiameter+1.75*G);

%% Rate Related Constants
gasConst = .008314; %kj/mol-K
k1 = 10^(4.2)*exp(-40.1/gasConst/T); %mol/(L-catalyst-hr)kPa^1.5
rxn1 = k1*F_C2H4*F_Cl2^0.5/F_tot^1.5*P^1.5;

k2 = 10^13.23*exp(-128.04/0.008314/T); 
rxn2 = k2*F_Prod*F_Cl2^0.5/F_tot^1.5*P^1.5;

k3 = 10^6.78*exp(-112/.008314/T);
rxn3=k3*F_O2*F_Cl2^0.5*F_C2H4/F_tot^2.5*P^2/5;

k4f = 1000*exp(17.13-13000/1.987/T);
k4b = exp(5.4+16000/1.987/T);
K4 = k4f/k4b;
rxn4 = K4*F_O2/F_Cl2;

%% The ODEs
% mass balance
rC2H4 = -rxn1-rxn3;
rProd = -rxn2 + rxn1;
rHCl = -2*rxn1 - rxn2 - 2*rxn4;
rO2 = -0.5*rxn1-0.5*rxn2-3*rxn3-0.5*rxn4;
rCO2 = 2*rxn3;
rH2O = rxn1+rxn2+2*rxn3+rxn4;
rCl3Eth = rxn2;
rCl2 = rxn4;

%Ergun equation
rP = -beta_0/(1-phi)/Ac*(P_0/P)*(T/T_0)*(F_tot/F_tot_0);
%% Convert to original script
dC2H4 = rC2H4;
dProd = rProd;
dHCl = rHCl;
dO2 = rO2;
dCO2 = rCO2;
dH2O = rH2O;
dCl3Eth = rCl3Eth;
dCl2 = rCl2;

dP = rP;

Homework8_ODE = [dProd; dC2H4; dHCl; dO2; dCO2; dH2O; dCl3Eth; dCl2; 0; 0; dP];
end 
