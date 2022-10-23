% This code is for CHEN 4520 Chemical Process Synthesis (Fall 2022)
% Homework 8
% Written by: Aimee Maravi, Abigail Hutabarat 
% Sean Oishi-Holder, and Mohammed Alahmed

clear all
clc
clear figure

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
Tube_temp_in= 466.32; %assuming the reaction happens at 533K
Coolant_temp_in = 300;
Pressure_i = 2000; %inlet pressure in kPa
IC = [F_Product_in, F_C2H4_i, F_HCl_i, F_O2_i, F_CO2_i, F_H2O_i, F_Cl3Eth_i, F_Cl2_i, Tube_temp_in, Coolant_temp_in, Pressure_i];
%% this is my Domain
V_I = 0;
%guessing the length is of the reactor is z? do we also just guess and
%check this value???????????????/
V_Domain = [V_I 6]; %Define the volume domain
%% SOLVE ODE 
[Vsol, Ysol] = ode45('DesignODETake3', V_Domain, IC);

%% Data handling 
% Extraction Solutions for each variable 
Product_sol = Ysol (:, 1) ;
C2H4_sol = Ysol(:,2); 
HCl_sol = Ysol(:,3);
O2_sol = Ysol(:,4);
CO2_sol = Ysol(:,5);
H2O_sol = Ysol(:,6);
Cl3Eth_sol = Ysol(:,7); 
Cl2_sol = Ysol(:,8);
TubeTemp_sol = Ysol(:,9);
CoolantTemp_sol = Ysol(:,10);
Pressure_sol = Ysol(:,11);

%% plotting 

fig1 = figure(1);
set(fig1, 'name', 'FlowRates');
%I realized we need a reactor volume rate
xlabel('Reactor Volume (m3)');
ylabel('FlowRates');
hold on
plot(Vsol,Product_sol, 'LineWidth', 4)
plot(Vsol,C2H4_sol, 'LineWidth', 4)
plot(Vsol,HCl_sol, 'LineWidth', 4)
plot(Vsol,CO2_sol, 'LineWidth', 4)
plot(Vsol,O2_sol, 'LineWidth', 4)
plot(Vsol,H2O_sol, 'LineWidth', 4)
plot(Vsol,Cl3Eth_sol, 'LineWidth', 4)
plot(Vsol,Cl2_sol, 'LineWidth', 4, 'Color','magenta')
hold off
legend('C2H4Cl2','C2H4','HCl','CO2','O2','H2O','C2H3Cl3','Cl2')

fig2 = figure(2);
set(fig2,'name', 'Reactor Temperature Profile')
area(Vsol, TubeTemp_sol);
xlabel('Tube Volume');
ylabel('Temperature ');

fig3 = figure(3); 
set(fig3,'name', 'mAb')
area(Vsol, CoolantTemp_sol);
xlabel('Reactor Volume');
ylabel('Coolant Temperature');

fig4 = figure(4);
set(fig4, 'name','pressure')
plot(Vsol,Pressure_sol);
xlabel('Reactor Volume (m3)');
ylabel('pressure (kPa)');

