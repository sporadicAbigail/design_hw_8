function [Homework8_ODE] = DesignODETake3(V,Y)

%% Reactor Volume Parameters
% this comes first because everything is going to be done on a per tube basis
numTubes = 1; %CHANGE ME
reactorLength = 4; %meters, CHANGE ME
mDotC = 30000; %kg/hr

%% Other things that change
T = Y(9); %Kelvin
T_0 = 500; %Kelvin, it works with 490
P = Y(11); %kPa
P_0 = 2000; %kPa

Tc = Y(10); %Kelvin

F_Prod = Y(1)/numTubes; %product
F_C2H4 = Y(2)/numTubes; %
F_HCl = Y(3)/numTubes; %
F_O2 = Y(4)/numTubes; %
F_CO2= Y(5)/numTubes;
F_H2O =Y(6)/numTubes;
F_Cl3Eth = Y(7)/numTubes;
F_Cl2 =Y(8)/numTubes;

F_tot = F_Prod + F_C2H4 + F_HCl + F_O2 + F_CO2 + F_H2O + F_Cl3Eth + F_Cl2; 

%% Rate Related Constants
gasConst = .008314; %kj/mol-K
k1 = 10^(4.2)*exp(-40.1/gasConst/T); %mol/(L-catalyst-hr)kPa^1.5
rxn1 = k1*F_C2H4*F_Cl2^0.5/F_tot^1.5*P^1.5;

k2 = 10^13.23*exp(-128.04/0.008314/T); 
rxn2 = k2*F_Prod*F_Cl2^0.5/F_tot^1.5*P^1.5;

k3 = 10^6.78*exp(-112/.008314/T);
rxn3=k3*F_O2*F_Cl2^0.5*F_C2H4/F_tot^2.5*P^2.5;

k4f = 1000*exp(17.13-13000/1.987/T);
k4b = exp(5.4+16000/1.987/T);
K4 = k4f/k4b;
rxn4 = K4*F_O2/F_Cl2;

%% Ergun Related Equations

F_tot_0_overall = 970240 + 2910750 + 970240 + 5820; %[mol/hr], constant
F_tot_0 = F_tot_0_overall/numTubes; %[mol/hr] per tube, constant

rho_0 = 1.002; %kg/m3, constant
rho = rho_0*(P/P_0)*(T_0/T)*(F_tot_0/F_tot); %kg/m3, changing

volumetricFlowRate_tot_0 = 166430.8*1000; %L/hr, calculated in writeup, overall, constant
volumetricFlowRate_0 = volumetricFlowRate_tot_0/numTubes; %L/hr, per tube, constant

Tau = 3.6; %space time in seconds, given 
reactorVol = Tau/3600*volumetricFlowRate_tot_0 %3600 to convert from s to hr: [L]
reactorVolm3 = reactorVol/1000; %[m3]
Ac = reactorVolm3/numTubes/reactorLength; %[m]
tubeDiameter=2*sqrt(Ac/pi) %cross sectional area [m2]
%here, Matlab calculates a reactor length
% intentionally without a semicolon so it comes out in our command window
% as a result


superFicVelocity = volumetricFlowRate_0/1000/Ac/3600; % m/s
G = rho*superFicVelocity; % superficial mass velocity [kg/m2-s]

phi = 0.50; %void fraction/porosity [unitless], constant
mu = 0.023*1000; %viscosity of gas [kPa-s] or [1000 kg/m-s], constant
gc = 1; %unit conversion [metric], constant
particleDiameter = tubeDiameter/8; %based on heuristic [m], constant

beta_0 = G*(1-phi)/(rho_0*gc*particleDiameter*phi^3)*(150*(1-phi)*mu/particleDiameter+1.75*G);

%% Energy Balance Related Equations
C_C2H4 = [0.3338*10^(5) 0.9479*10^5 1.596*10^3 0.551*10^5 740.8]; %array of constants for C2H4
Cp_C2H4 = C_C2H4(1) + C_C2H4(2)*(C_C2H4(3)/T/sinh(C_C2H4(3)/T))^2 +C_C2H4(4)*(C_C2H4(5)/T/cosh(C_C2H4(5)/T))^2; %kj/mol-K
C_HCl = [0.29157*10^5 0.09048*10^5 2.0938*10^3 -0.00107*10^5 120]; %array of constants for HCl
Cp_HCl = C_HCl(1) + C_HCl(2)*(C_HCl(3)/T/sinh(C_HCl(3)/T))^2 +C_HCl(4)*(C_HCl(5)/T/cosh(C_HCl(5)/T))^2; %kj/mol-K
C_O2 = [0.29103*10^5 0.1004*10^5 2.5265*10^3 0.09356*10^5 1153.8]; %array of constants for O2
Cp_O2 = C_O2(1) + C_O2(2)*(C_O2(3)/T/sinh(C_O2(3)/T))^2 +C_O2(4)*(C_O2(5)/T/cosh(C_O2(5)/T))^2; %kj/mol-K
C_prod = [0.5521*10^5 1.205*10^5 1.502*10^3 0.8719*10^5 653.5]; %array of constants for prod
Cp_Prod = C_prod(1) + C_prod(2)*(C_prod(3)/T/sinh(C_prod(3)/T))^2 +C_prod(4)*(C_prod(5)/T/cosh(C_prod(5)/T))^2; %kj/mol-K
C_H2O = [0.33363*10^5 0.2679*10^5 2.6105*10^3 0.08896*10^5 1169]; %array of constants for H2O
Cp_H2O = C_H2O(1) + C_H2O(2)*(C_H2O(3)/T/sinh(C_H2O(3)/T))^2 +C_H2O(4)*(C_H2O(5)/T/cosh(C_H2O(5)/T))^2;
C_Cl3Eth = [0.66554*10^5 1.1257*10^5 1.5454*10^3 0.97196*10^5 717.04]; %array of constants for Cl3Eth
Cp_Cl3Eth = C_Cl3Eth(1) + C_Cl3Eth(2)*(C_Cl3Eth(3)/T/sinh(C_Cl3Eth(3)/T))^2 +C_Cl3Eth(4)*(C_Cl3Eth(5)/T/cosh(C_Cl3Eth(5)/T))^2; 
C_CO2 = [0.2937*10^5 0.3453*10^5 1.428*10^3 0.264*10^5 588]; %array of constants for CO2
Cp_CO2 = C_CO2(1) + C_CO2(2)*(C_CO2(3)/T/sinh(C_CO2(3)/T))^2 +C_CO2(4)*(C_CO2(5)/T/cosh(C_CO2(5)/T))^2; %kj/mol-K
C_Cl2 = [0.29142*10^5 0.09176*10^5 0.949*10^3 0.1003*10^5 425]; %array of constants for Cl2
Cp_Cl2 = C_Cl2(1) + C_Cl2(2)*(C_Cl2(3)/T/sinh(C_Cl2(3)/T))^2 +C_Cl2(4)*(C_Cl2(5)/T/cosh(C_Cl2(5)/T))^2; %kj/mol-K
Cp_dowtherm = 1.125 + 0.0025*(T); %kj/kg-K

%should these be neg or pos
hrxn1 = -239.111; %kj/mol
hrxn2 = -162.091; %kj/mol
hrxn3 = -1323.155; %kj/mol
hrxn4 = -228.786; %kJ/mol

Ua = 300*tubeDiameter*pi*reactorLength/(reactorVolm3/numTubes)*0.0036; % [kJ/(hr L-cat K) Heat capacity*surface area of heat transfer/volume

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

%Thermal equation
numerator = (-rxn1*hrxn1-rxn2*hrxn2-rxn3*hrxn3-rxn4*hrxn4)-Ua*(T-Tc);
denominator = ((1-phi)*(F_C2H4*Cp_C2H4+F_Prod*Cp_Prod+F_HCl*Cp_HCl+F_O2*Cp_O2+F_CO2*Cp_CO2+F_H2O*Cp_H2O+F_Cl3Eth*Cp_Cl3Eth+F_Cl2*Cp_Cl2));
rT = numerator/denominator;

%cooling thermal balance
rTc = Ua*(T-Tc)/(mDotC/numTubes)/Cp_dowtherm;

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

dT = rT;

dTc = rTc;
Homework8_ODE = [dProd; dC2H4; dHCl; dO2; dCO2; dH2O; dCl3Eth; dCl2; dT; dTc; dP];
% %Goal
% F_C2H4_fin = 0.04*970240/numTubes;
% if F_C2H4 <= F_C2H4_fin 
%     return;
% end
end 
