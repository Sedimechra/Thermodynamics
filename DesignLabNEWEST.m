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

Temp = 303; % in K
TempSpace = 0; 

%((Arad * EmissivityRadiator * IRSolarUse) + IRSpaceCraft + (Arad * EmissivityRadiator * IRSpaceCraft)) = (Arad * EmissivityRadiator)

AradSummer = -20/((IRSolarUseSummer*AbsorbanceRadiator*cosd(phi)) + (IRSpaceCraftSummer * EmissivityRadiator) - (EmissivityRadiator *(5.67*10^-8)* (Temp ^4 - TempSpace ^4)))
SideRad = sqrt(AradSummer)


AradWinter = -20/((IRSolarUseWinter*AbsorbanceRadiator*cosd(phi)) + (IRSpaceCraftWinter * EmissivityRadiator) - (EmissivityRadiator *(5.67*10^-8)* (Temp ^4 - TempSpace ^4)))
SideRad = sqrt(AradWinter)

IRSpaceCraftEquinox = (63+88)/2; %W/m^2
IRSpaceCraftEclipse = 11;
Aradeq = -20/((IRSolar*AbsorbanceRadiator) + (IRSpaceCraftEquinox * EmissivityRadiator) - (EmissivityRadiator *(5.67*10^-8)* (Temp ^4 - TempSpace ^4)))
%fprintf("IR Use is %f", IRSolarUse)


SolarFlux = [];
TempSpaceCraft = [];
TempSpaceCraft2 = [];
OpPower = [];
EPower = [];
for i = 0:0.01:2*pi
    if i <= pi 
       
        SolarF = Aradeq * IRSolarUseWinter * AbsorbanceRadiator *cosd(23.5)* sin(i); %+ (IRSpaceCraftWinter*AradWinter * AbsorbanceRadiator) ;
    
        OpP = 74.0071 - SolarF;
        if OpP < 0
           OpP = 0; 
        end
        
    else 
        SolarF = 0;
    end 
    SolarFlux = [SolarFlux;SolarF];
    OpPower = [OpPower; OpP];
    maxind1 = max(SolarFlux);
    T = ((SolarF + 20 + (Aradeq * EmissivityRadiator*IRSpaceCraftWinter))/(Aradeq * EmissivityRadiator *(5.67*10^-8)))^(1/4);
    TempSpaceCraft = [TempSpaceCraft;T];
    T2 = ((SolarF + OpP + 20 + (Aradeq * EmissivityRadiator*IRSpaceCraftWinter))/(Aradeq * EmissivityRadiator *(5.67*10^-8)))^(1/4);
    TempSpaceCraft2 = [TempSpaceCraft2;T2]; 
end 

%emergency power for winter can just be normal 20
for j = 1:629
     EP = 20;
     EPower = [EPower; EP];
end


figure  
plot(linspace(0,24,629), SolarFlux, 'LineWidth', 3)
xlabel("Time (Hours, beginning at noon)")
ylabel("Solar Flux [W]")
title('Solar Flux During Winter Solstice')
grid on
grid minor


figure
yyaxis left
plot(linspace(0,24,629), TempSpaceCraft, 'b--', 'LineWidth', 3)
xlabel("Time (Hours, beginning at noon)")
ylabel("Temperature [K]")
title('Unheated Temperature of the Radiator During Winter Solstice')
grid on
grid minor
hold on
plot(linspace(0,24,629), TempSpaceCraft2, 'b', 'LineWidth', 3)
yyaxis right
ylabel('Operational Power [W]')
plot(linspace(0,24,629), OpPower, 'r', 'LineWidth', 3)
legend('Unheated Temperature', 'Heated Temperature')



figure
yyaxis left
plot(linspace(0,24,629), TempSpaceCraft, 'b', 'LineWidth', 3)
xlabel("Time (Hours, beginning at noon)")
ylabel("Temperature [K]")
title(' Temperature of the Radiator During Winter Solstice')
grid on
grid minor
hold on
plot(linspace(0,24,629), TempSpaceCraft2, 'b', 'LineWidth', 3)
yyaxis right
ylabel('Operational Power [W]')
plot(linspace(0,24,629), OpPower, 'r', 'LineWidth', 3)
legend('Unheated Temperature', 'Heated Temperature')

figure
plot(linspace(0,24,629), OpPower, 'k', 'LineWidth', 3)
xlabel("Time (Hours, beginning at noon)")
ylabel("Power [W]")
title('Power During Winter Solstice')
grid on
grid minor
hold on
plot(linspace(0,24,629), EPower, 'm', 'LineWidth', 3)
legend('Operational Power','Emergency Power')


SolarFlux = [];
TempSpaceCraft = [];
OpPower = [];
EPower = [];
for i = 0:0.01:2*pi
    if i <= pi 
        SolarF = Aradeq * IRSolarUseSummer * AbsorbanceRadiator * cosd(23.5)*sin(i); %+ (IRSpaceCraftWinter*AradWinter * AbsorbanceRadiator) ;
         OpP = 67 - SolarF;
        if OpP < 0
           OpP = 0; 
        end
        
    
    else 
        SolarF = 0;
    end 
    SolarFlux = [SolarFlux;SolarF];
    OpPower = [OpPower; OpP];
    maxind2 = max(SolarFlux);
    T = ((SolarF + 20 + (Aradeq * EmissivityRadiator*IRSpaceCraftSummer))/(Aradeq * EmissivityRadiator *(5.67*10^-8)))^(1/4);
    TempSpaceCraft = [TempSpaceCraft;T];
    
end 

%emergency power calculations
for j = 1:629
     if TempSpaceCraft(j) > (273-40)
        EP = 20;
     else
         EP = 30;
     end
     EPower = [EPower; EP];
end
     
     

figure
plot(linspace(0,24,629), TempSpaceCraft, 'r', 'LineWidth', 3)
xlabel("Time (Hours, beginning at noon)")
ylabel("Temperature [K]")
title('Unheated Temperature of the Radiator During Summer Solstice')
grid on
grid minor

figure  
plot(linspace(0,24,629), SolarFlux, 'r', 'LineWidth', 3)
xlabel("Time (Hours, beginning at noon)")
ylabel("Solar Flux [W]")
title('Solar Flux During Summer Solstice')
grid on
grid minor

figure
plot(linspace(0,24,629), OpPower, 'k', 'LineWidth', 3)
xlabel("Time (Hours, beginning at noon)")
ylabel("Power [W]")
title('Power During Summer Solstice')
grid on
grid minor
hold on
plot(linspace(0,24,629), EPower, 'm', 'LineWidth', 3)
legend('Operational Power','Emergency Power')






% 
SolarFlux = [];
TempSpaceCraft = [];
OpPower = [];
EPower = [];
for i = 0:0.01:2*pi
    if i <= 2.98 
        SolarF = Aradeq * IRSolar * AbsorbanceRadiator * sin(i); %+ (IRSpaceCraftWinter*AradWinter * AbsorbanceRadiator) ;
        T = ((SolarF + 20 + (Aradeq * EmissivityRadiator*IRSpaceCraftEquinox))/(Aradeq * EmissivityRadiator *(5.67*10^-8)))^(1/4);
        OpP = 76 - SolarF;
        if OpP < 0
           OpP = 0; 
        end
        
    elseif i > 2.98 && i < 3.30
        
        SolarF = 0;
        
        T = ((SolarF+(Aradeq * EmissivityRadiator*IRSpaceCraftEclipse))/(Aradeq * EmissivityRadiator * (5.67*10^-8)))^(1/4);
       OpP = 76 - SolarF;
        if OpP < 0
           OpP = 0; 
        end
    
    else
       SolarF = 0;
        T = ((SolarF + 20 + (Aradeq * EmissivityRadiator*IRSpaceCraftEquinox))/(Aradeq * EmissivityRadiator *(5.67*10^-8)))^(1/4);
        OpP = 76 - SolarF;
        if OpP < 0
           OpP = 0; 
        end
    end
    SolarFlux = [SolarFlux;SolarF];
    OpPower = [OpPower; OpP];
   % T = ((SolarF + 20 + (Aradeq * EmissivityRadiator*IRSpaceCraftEquinox))/(Aradeq * EmissivityRadiator *(5.67*10^-8)))^(1/4);
    TempSpaceCraft = [TempSpaceCraft;T];
    maxind3 = max(SolarFlux);
end 


%emergency power calculations
for j = 1:629
     if TempSpaceCraft(j) > (273-40)
        EP = 20;
     elseif TempSpaceCraft(j) < (273-50)
         EP = 50; %gotta change this to account for super low temp
     else
         EP = 30;
     end
     EPower = [EPower; EP];
end


figure  
plot(linspace(0,24,629), SolarFlux, 'g', 'LineWidth', 3)
xlabel("Time(Hours, beginning at noon)")
ylabel("Solar Flux [W]")
title('Solar Flux During Equinox')
grid on

figure
plot(linspace(0,24,629), TempSpaceCraft, 'g', 'LineWidth', 3)
xlabel("Time(Hours, beginning at noon)")
ylabel("Temperature [K]")
title('Unheated Temperature of the Radiator During Equinox')
grid on

figure
plot(linspace(0,24,629), OpPower, 'k', 'LineWidth', 3)
xlabel("Time (Hours, beginning at noon)")
ylabel("Power [W]")
title('Power During Equinox')
grid on
grid minor
hold on
plot(linspace(0,24,629), EPower, 'm', 'LineWidth', 3)
legend('Operational Power','Emergency Power')

