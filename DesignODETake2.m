function [Homework8_ODE] = DesignODETake2(V,Y)
%% output variables 

%Chemical Species molar flow rates
F_Cl2Eth = Y(1); 
F_C2H4 = Y(2); %
F_HCl = Y(3); %
F_O2 = Y(4); %
F_CO2= Y(5);
F_H2O =Y(6);
F_Cl3Eth = Y(7);
F_Cl2 =Y(8);

%Changing F total
F = F_Cl2Eth+F_C2H4+F_HCl+F_O2+F_CO2+F_H2O+F_Cl3Eth+F_Cl2; 

%temperature of the tube and the coolant
T = Y(9);
Tc = Y(10);
rT = 0; %change in temperature

%kPa go back and check this
P= Y(11); %pressure in 
rP = 0; % just in  case fam
%Homework8_ODE = [d_product; d_C2H4; d_HCl; d_O2; d_CO2; d_H2O; d_Cl3Eth; d_Cl2; d_T_Tube; d_T_Coolant;d_P];
%% Defining Constants 

%THE GUESSES
% The following values are values that we choose/guess
% reactor length is also unknown, but we calculate
numTubes = 1; %this must be an integer
tubeDiameter = 1; %We will decide this here, because we have some 
% engineering intuition here

%this is initial total volumetric flow rates combining all the tubes [m^3/hr]
Product_in_vol = 0; 
C2H4_i_vol= 1729.38;
HCl_i_vol = 160943.2;
O2_i_vol = 3724.53;
CO2_i_vol = 0;
H2O_i_vol = 0;
Cl3Eth_i_vol = 0;
Cl2_i_vol = 33.789;
%Inital volumetric flow rate in
VolFlow_i = C2H4_i_vol + HCl_i_vol + O2_i_vol + CO2_i_vol + H2O_i_vol + Cl3Eth_i_vol + Cl2_i_vol; 

%this is initial molar flow rates combining all the tubes [mol/hr]
%We are setting these 
F_Product_in = 0;
F_C2H4_i= 970240; %[mol/hr]
F_HCl_i = 2910750;
F_O2_i= 970240;
F_CO2_i = 0;
F_H2O_i = 0;
F_Cl3Eth_i = 0;
F_Cl2_i = 5820;
%initial molar flow rate  is constant
molFlow_i = F_C2H4_i + F_HCl_i + F_O2_i + F_CO2_i + F_H2O_i + F_Cl3Eth_i + F_Cl2_i; 
% divide by the number of tubes so it's a single tube
F_0_single = molFlow_i/numTubes;

Tau = 3.6; %space time in seconds 
reactorVol = Tau/3600*VolFlow_i; %3600 to convert from s to hr
%cross sectional area
Ac=pi*(tubeDiameter/2)^2;
%here, Matlab calculates a reactor length
% intentionally without a semicolon so it comes out in our command window
% as a result
reactorLength = reactorVol/numTubes/Ac;


%heat capacities 
CP_C2H4=1;
CP_Cl2Eth=1;
CP_HCl=1;
CP_O2=1;
CP_CO2=1;
CP_H2O=1;
CP_Cl3Eth=1;
CP_C12=1;
CP_C= 2.25; %Coolant heat capacity, dowtherm A

%other constants guess and check 
%V= 1; %Volumetric velocity
%G = 1 ;%gravity Not Constant
g_c = 1;%some unit conversion (in metric units)
D_p = tubeDiameter/8;%bead diameter, based on heuristic
Void = 0.50; %void fraction
mu = 0.023; %gas viscosity, assume constant throughout the pipe [mPa s]
M_DotC = 1; %mass flow of the coolant
U = 300; %Overall heat transfer coefficient [W/(m^2-K)]
A = tubeDiameter*pi*reactorLength/reactorVol; %surface area of heat exchange/volume reactor
rho_0 = 1.002; %initial gas density for the ergun equation NOT constant [kg/m3]
rho_bulk = 1975; %kg/m3
rho_cat = 3950; %kg catalyst/m3
R_gas = 1; % CHECK THE UNITS ON R GAS 
P_inlet = 2000; % kPa 

%% Heats of Reaction [kJ/mol]
H_1 = -239.111;
H_2 = -162.090;
H_3 = -1323.155;
H_4 = 228.8;

%% Base rate equations 
K_1 = 10^(4.2)*exp(-40.1/(8.314*T)); %K_1 temperature Arrhenius equation
R_1 = K_1*(F_C2H4*(F_Cl2^0.5)/(F^(1.5)))*P^(1.5);

K_2 = 10^(13.23)*exp(-128.04/(8.314*T));
R_2 = K_2*(F_Cl2Eth*(F_Cl2^0.5)/(F^1.5))*P^1.5;

K_3 =10^(6.78)*exp(-112/(8.314*T));
R_3 = K_3*(F_O2*(F_Cl2^0.5)*(F_C2H4)/(F^2.5))*P^2.5;

K_4f = 1000*exp(5.4+(160000/(1.987*T)));
K_4b = exp(5.4+(160000/(1.987*T)));
%*****we are using the experimental values found in the paper not the 2
%termed rate equation
K_4 = K_4f/K_4b;
R_4 = K_4*((F_O2)/(F_Cl2)); %check the units for k4 and make sure that works out

%% Species balances
rC2H4 = -(R_1 + R_3);
rCl2Eth = -(R_2-R_1);
rHCl = -(2*R_2 + R_1 + 2*R_4);
rO2 = -(0.5*R_1+ 0.5*R_2+ 3*R_3 + 0.5*R_4);
rCO2 = -(-2*R_3);
rH2O = -(-R_1 -R_2 -2*R_3 -R_4);
rCl3Eth = R_2;
rCl2= R_4;
rF = rC2H4+ rCl2Eth + rHCl + rO2 + rCO2 + rH2O + rCl3Eth + rCl2; %total flow rate based off of species balance (changes)

%% Thermal Balances 
%per tube 

%Coolant fluid

%% Pressure Considerations %which we also have to guess and check.....

% FIX THIS P

%% the end?
d_product = rCl2Eth; % [mol/Lcat-hr]
d_C2H4 = rC2H4;
d_HCl = rHCl;
d_O2 = rO2;
d_CO2 = rCO2;
d_H2O = rH2O;
d_Cl3Eth = rCl3Eth;
d_Cl2 = rCl2;
d_T_Tube = rT;
d_T_Coolant = 0;
d_P = rP;

Homework8_ODE = [d_product; d_C2H4; d_HCl; d_O2; d_CO2; d_H2O; d_Cl3Eth; d_Cl2; d_T_Tube; d_T_Coolant;d_P];
end 