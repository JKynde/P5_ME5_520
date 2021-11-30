function [param] = StructCreator()
% This script creates a struct will all relevant parameters, such that it can be passed to a function.
%% Editable running options
% Kan ændres efter structen er passet til det pågældende scribt som skal
% bruge structen.
param = struct;
param.Tlmode=1;
param.mode = 3;
%% Non-calculated parameters
% Priezoelektriske Materiale konstanter (PIC181).
param.rho_P = 7800; % Densitet i [kg/m^3]
param.epsilon_0 = 8.854*10^(-12); % permittivity of free space [F/m]
param.PermittivityPolarizationInDirection = 1200; % Relative permittivity in the polarization direction = epsilon_33^T/epsilon_0
%PermittivityNormalDirection = 1500; % Relative permittivity in the norm direction epsilon_11^T/epsilon_0


% Piezoelectric Charge Coefficients [10^-12C/N)
param.d_33 = 265*10^(-12);
    % Acousto-mechenical properties
%Elastic compliance coefficients [10^-12m^2/N]
param.S33E = 14.2*10^(-12);
param.c33D = 16.6*10^10; % Elastic stiffness coefficient [10^10N/m^2]
param.C_0 = 0.104*10^(-9); % Kapacitans af pladen i farad 
%c33D = ( d_33 / ( S33E * epsilon33T) )^2*epsilon33T + 1 / (S33E);
%Alternativ måde at regne c33D på

% Der er en antagelse at h_33 i Ultrasonic Nondestructive Evaluation
% systems er den samme som i den artiklen hvor formlen for h_33 er fundet.

%Geometri af en enkelt piezoplade
param.d_p = 0.002; %Plate thickness [m]
param.r_transducer = 0.005/2;
param.d_heads = 44.8*10^-3; % Distance between transducer heads.

%Backing
param.ElasticModolusBacking = 2.1*10^11; %Backing materiale konstanter.[N/m^2]
param.rho_b = 7800;%Backing densitet [kg/M^3]m

%Front layer parameters. Altså housingen ude foran.
param.l_a = 3.2*10^(-3); %layer thickness [m]
param.rho_f = 2700; %Densitet af front layer. [kg/m^3]
param.ElasticModolusFront = 6.9*10^10; %Young's modolus for front end materialet

%adhesive parameters
param.l_aa = 0.015*10^(-3);
param.rho_a = 1146;
param.ElasticModolusAdhesive = 5.9*10^10;

% %Mellem layer parameters. Er kobber nu, og ikke plain carbon steel.
param.l_m =0.5*10^(-3);
param.rho_m = 8933;
param.ElasticModolusMellem = 1.21*10^11;

%Partikle parameters:
%luft
% param.rho_p = 1.225; %kg/m^3 %rho luft = 1.225
% param.v_0p = 340.19744;% m/s eller 761 i retard units 
% brug disse for polysterene(tunge partikler)
param.rho_p = 1050; 
param.v_0p = 1700; % og denne

% Oil/fuidparameters parameters ISO VG-32 ISO.
% param.rho_oil = 857; % Densitet af olien [kg/m^3]
% param.v_0Oil = 1493; % Lydens hastighed i olien.
% param.nu_oil = 32*10^-6; % Oliens kinematiske viskøsitet  https://eurol.com/da/produkter/eurol-synmax-pao-iso-vg-320/
%Water
param.rho_oil = 1000; % Densitet af vand [kg/m^3]
param.v_0Oil = 1480; % Lydens hastighed i Vand.
param.nu_oil = 1*10^-6; % Vand kinematiske viskøsitet  https://eurol.com/da/produkter/eurol-synmax-pao-iso-vg-320/


%Cabel parameters
 param.l_c = 0.85; %Længde af kabel
 param.b_c = 0.5*1.5*10^(-3); %Ydre radius af kabel
 param.ElasticModolusCabel = 11.7*10^10; % Elastic Modolus af kobber
 param.rho_c = 8960; %kg/m3 af kobber
 param.sigma_c = 58.5*10^6; %conductivity of copper
 param.mu_c = 1.256629*10^-6; %permeability of copper
 param.c_c = 2*10^8; %speed of signal in copper
 
%% Calculated parameters
param.Area = (param.r_transducer)^2*pi; % Plate area
param.epsilon33T = param.PermittivityPolarizationInDirection*param.epsilon_0; %epsilon 33T skal bruges senere til at udregnes h_33, som vi skal bruge.
param.h_33=param.d_33/(param.S33E*param.epsilon33T); % Formel for h_33, som vi tror er Piezoelectric stiffness constant er regnet ud fra side 15 i https://link.springer.com/content/pdf/bbm%3A978-94-007-0579-1%2F1.pdf

param.lambda = param.v_0Oil/10^6; % Bølgelængde v/f ved 10 MHz

param.v_0 = sqrt(param.c33D/param.rho_P); % wave speed of compressional waves in the piezoelectric plate
param.v_0m = sqrt(param.ElasticModolusMellem/param.rho_m); % Wave speed of compressional waves in the mellemlag.. nice danlish
param.v_0b = sqrt(param.ElasticModolusBacking/param.rho_b); % Wave speed of compressional waves in the backing material
param.v_0f = sqrt(param.ElasticModolusFront/param.rho_f); % Wave speed of compressional waves in the front layer
param.v_0a = sqrt(param.ElasticModolusAdhesive/param.rho_a); % Wave speed of compressional waves in the adhesive
param.v_0c = sqrt(param.ElasticModolusCabel/param.rho_c); % Wave speed of compressional waves in the adhesive

param.n=param.h_33*param.C_0; % ElectricMatrix parameter.

param.Z0a = param.rho_P*param.v_0*param.Area; % plane wave acoustic impedance of the piezoelectric plate.
param.Zba = param.rho_b*param.v_0b*param.Area; % <- vides ikke endnu, men det er: Corresponding acoustic impedance of the backing Z_b^A Vi gætter foreløbigt på =rho_backing*v_0,backing*S Problemet er at der står den er afhængig af frekvens. Plus den er cone-shaped.
param.Zma = param.rho_m*param.v_0m*param.Area; % Acoustic impedance mellemlag
param.Zaa = param.rho_a*param.v_0a*param.Area; % Acoustic impedance adhesive.
param.Z0a_f = param.rho_f*param.v_0f*param.Area; % Acoustic impedance front layer ;
param.ZrAa =param.v_0Oil*param.rho_oil*param.Area; % Acoustic impedance of radiating medium. In this case the oil

param.kappa_l=1/(param.rho_oil*param.v_0Oil^2); %Gorkov lort
param.kappa_p = 1/(param.rho_p*param.v_0p^2); % mere gorkov lort

param.R_reflect = (param.rho_f*param.v_0f-param.rho_oil*param.v_0Oil) / (param.rho_oil*param.v_0Oil + param.rho_f*param.v_0f); % Reflection coefficient between transducer face and oil.

param.f_1 = 1 - param.kappa_p/param.kappa_l; % En konstant der skal bruges i gorkov
param.f_2 = (2*(param.rho_p-param.rho_oil)) / (2*param.rho_p + param.rho_oil); % Også en konstant der skal bruges i gorkov

param.Contrast_factor = (5*param.rho_p-2*param.rho_oil) / (2*param.rho_p+param.rho_oil) - param.kappa_p/param.kappa_l; % Siger noget om hvor vidt partiklerne samles i nodes eller antinodes

param.mu_oil = param.nu_oil*param.rho_oil; % nu fluid

end

