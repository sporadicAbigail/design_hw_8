% This code is for CHEN 4520 Chemical Process Synthesis (Fall 2022)
% Homework 8
% Written by: Aimee Maravi, Abigail Hutabarat, Mady Murphy, 
% Sean Oishi-Holder, and Mohammed Alahmed

clear
clc

T_out = 533; %outlet Temperature in kelvin 
Cbd = 1975; %kg cat/m^3 this is catalyst bulk density 
Void = 0.50; %catalyst voidage
%heurisitc to account for later: catalyst particle diameter is an eigth of
%the tube diameter

%% Initial Conditions
%T0 = 1; %what units 

%this is molar flow rates [mol/hr]
F_Product_in = 0;
F_C2H4_i= 970240;
F_HCl_i = 2910750;
F_O2_i= 970240;
F_CO2_i = 0;
F_H2O_i = 0;
F_Cl3Eth_i = 0;
F_Cl2_i = 5820;
%initial molar flow rate  is constant
molFlow_i = F_C2H4_i + F_HCl_i + F_O2_i + F_CO2_i + F_H2O_i + F_Cl3Eth_i + F_Cl2_i; 

%Guess T initial
Tube_temp_in=0;
Coolant_temp_in = 0;
Pressure_i = 2000; %inlet pressure in kPa
IC = [F_Product_in, F_C2H4_i, F_HCl_i, F_O2_i, F_CO2_i, F_H2O_i, F_Cl3Eth_i, F_Cl2_i, Tube_temp_in, Coolant_temp_in, Pressure_i];
%Homework8_ODE = [d_product; d_C2H4; d_HCl; d_O2; d_CO2; d_H2O; d_Cl3Eth; d_Cl2; d_T_Tube; d_T_Coolant;d_P];
%% this is my Domain
V_I = 0;
%guessing the length is of the reactor is z? do we also just guess and
%check this value???????????????/
V_Domain = [V_I 1000]; %Define the temperature domain
%% SOLVE ODE 
[Vsol, Ysol] = ode45('DesignODE', V_Domain, IC);

%% Data handling 
% Extraction Solutions for each variable 
Product_sol = Ysol (:, 1) ;
C2H4_sol = Ysol(:,2); 
HCl_sol = Ysol(:,3);
O2_sol = Ysol(:,4);
CO2_sol = Ysol(:,5);
H2O_sol = Ysol(:,6);
Cl3Eth_sol = Ysol(:,7); 
F_Cl2_i = Ysol(:,8);
TubeTemp_sol = Ysol(:,9);
CoolantTemp_sol = Ysol(:,10);
Pressure_sol = Ysol(:,11);

%% plotting 

fig1 = figure(1);
set(fig1, 'name', 'FlowRates');
%I realized we need a reactor volume rate
area(Product_sol,Vsol);
xlabel('Reactor Volume (L)');
ylabel('FlowRates');
hold on 
plot(C2H4_sol)
hold on 
plot(HCl_sol)
hold on 
plot(O2_sol)
hold on 
plot(H2O_sol)
hold on
plot(Cl3Eth_sol)
hold on
plot(F_Cl2_i)

fig2 = figure(2);
set(fig2,'name', 'Reactor Temperature Profile')
area(Vsol, TubeTemp_sol);
xlabel('Tube Volume');
ylabel('Temperature ');

fig3 = figure(3); 
set(fig2,'name', 'mAb')
area(Vsol, CoolantTemp_sol);
xlabel('Reactor Volume');
ylabel('Coolant Temperature');
%we need reactor pressure- thats probably what one of the ODEs we're
%missing currentlly 

%
%points of confusion for me 1) what the heck is the units of Ysol?
%Ysol is solution vector for all the ODEs
%Xsol will be in units of volume, everything in terms of volume
%are we integrating by 2 variables? No - technically more than one