%% Design Lab 
clear; clc; close all; 

IRSolar = 1361; %at 1 AU in w/m^2
IRSpaceCraftWinter = 88;  %w/m^2
IRSpaceCraftSummer = 63;  %w/m^2
InstrumentHeat = 20; %w
Rmin = 0.9827; %in AU
Rmax = 1.0205; %in AU 
theta = 0; % orbit angle 
phi = 23.5; % inclination angle 

% EmissivityRadiator = 0.91; % +- 0.02 
% AbsorbanceRadiator = 0.15; % +- 0.02

%These are to test that our code is giving the correct output 
EmissivityRadiator = 0.85; % +- 0.02 
AbsorbanceRadiator = 0.2; % +- 0.02

R = Rmin; 
IRSolarUseSummer = 1361 * (1/Rmax)^2; 
IRSolarUseWinter = 1361 * (1/Rmin)^2; 

%IRTotal = IRSolarUse + IRSpaceCraft + InstrumentHeat; 

Temp = 293; % in K
TempSpace = 0; 

%((Arad * EmissivityRadiator * IRSolarUse) + IRSpaceCraft + (Arad * EmissivityRadiator * IRSpaceCraft)) = (Arad * EmissivityRadiator)

AradSummer = -20/((IRSolarUseSummer*AbsorbanceRadiator*cosd(phi)) + (IRSpaceCraftSummer * EmissivityRadiator) - (EmissivityRadiator *(5.67*10^-8)* (Temp ^4 - TempSpace ^4)))
SideRad = sqrt(AradSummer)


AradWinter = -20/((IRSolarUseWinter*AbsorbanceRadiator*cosd(phi)) + (IRSpaceCraftWinter * EmissivityRadiator) - (EmissivityRadiator *(5.67*10^-8)* (Temp ^4 - TempSpace ^4)))
SideRad = sqrt(AradWinter)

%fprintf("IR Use is %f", IRSolarUse)


SolarFlux = [];
TempSpaceCraft = [];
for i = 0:0.01:2*pi
    if i <= pi 
        SolarF = AradWinter * IRSolarUseWinter * AbsorbanceRadiator * sin(i); %+ (IRSpaceCraftWinter*AradWinter * AbsorbanceRadiator) ;
    else 
        SolarF = 0;
    end 
    SolarFlux = [SolarFlux;SolarF];
    T = ((SolarF + 20 + (AradWinter * EmissivityRadiator*IRSpaceCraftWinter))/(AradWinter * EmissivityRadiator *(5.67*10^-8)))^(1/4);
    TempSpaceCraft = [TempSpaceCraft;T];
end 

figure  
plot([1:629], SolarFlux)
xlabel("Time")
ylabel("Solar Flux [W]")

figure
plot([1:629], TempSpaceCraft)
xlabel("Time")
ylabel("Temperature [k]")



