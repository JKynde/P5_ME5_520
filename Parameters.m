%Definitions
j = sqrt(-1);

% Priezoelektriske Materiale konstanter (PIC181).
rho_P = 7800; % Densitet i [kg/m^3]
epsilon_0 = 8.854*10^(-12); % permittivity of free space [F/m]
PermittivityPolarizationInDirection = 1200; % Relative permittivity in the polarization direction = epsilon_33^T/epsilon_0
epsilon33T = PermittivityPolarizationInDirection*epsilon_0; %epsilon 33T skal bruges senere til at udregnes h_33, som vi skal bruge.
PermittivityNormalDirection = 1500; % Relative permittivity in the norm direction epsilon_11^T/epsilon_0


% Piezoelectric Charge Coefficients [10^-12C/N)
d_33 = 265*10^(-12);
    % Acousto-mechenical properties
%Elastic compliance coefficients [10^-12m^2/N]
S33E = 14.2*10^(-12);
c33D = 16.6*10^10; % Elastic stiffness coefficient [10^10N/m^2]
%c33D = ( d_33 / ( S33E * epsilon33T) )^2*epsilon33T + 1 / (S33E); %
%Alternativ måde at regne c33D på
C_0 = 0.104*10^(-9); % Kapacitans af pladen i farad 
h_33=d_33/(S33E*epsilon33T); % Formel for h_33, som vi tror er Piezoelectric stiffness constant er regnet ud fra side 15 i https://link.springer.com/content/pdf/bbm%3A978-94-007-0579-1%2F1.pdf
% Der er en antagelse at h_33 i Ultrasonic Nondestructive Evaluation
% systems er den samme som i den artiklen hvor formlen for h_33 er fundet.

%Geometri af en enkelt piezoplade
d_p = 0.002; %Plate thickness [m]
r_transducer = 0.005/2;
Area = (r_transducer)^2*pi; % Plate area
d_heads = 44.8*10^-3; % Distance between transducer heads.

%Backing
ElasticModolusBacking = 2.1*10^11; %Backing materiale konstanter.[N/m^2]
rho_b = 7800;%Backing densitet [kg/M^3]m

%Front layer parameters. Altså housingen ude foran.
l_a = 3.2*10^(-3); %layer thickness [m]
rho_f = 2700; %Densitet af front layer. [kg/m^3]
ElasticModolusFront = 6.9*10^10; %Young's modolus for front end materialet

%adhesive parameters
l_aa = 0.015*10^(-3);
rho_a = 1146;
ElasticModolusAdhesive = 5.9*10^10;

% %Mellem layer parameters.
l_m =0.5*10^(-3);
rho_m = 7800;
ElasticModolusMellem = 2.1*10^11;

%Partikle parameters:
%luft
rho_p = 1.225; %kg/m^3 %rho luft = 1.225
v_0p = 340.19744;% m/s eller 761 i retard units 
% brug disse for polysterene(tunge partikler)
% rho_p = 1050; 
% v_0p = 1700; % og denne

% Oil parameters ISO VG-32 ISO.
rho_oil = 857; % Densitet af olien [kg/m^3]
v_0Oil = 1493; % Lydens hastighed i olien.
nu_oil = 32*10^-6; % Oliens kinematiske viskøsitet  https://eurol.com/da/produkter/eurol-synmax-pao-iso-vg-320/
mu_oil = nu_oil*rho_oil;

%Cabel parameters
 l_c = 0; %Længde af kabel
 epsilon_c = 0; %Permativitet af kabel
 myh_c = 0; %Permabiliteten af kabel
 a_c = 0; %Indre radius af kabel
 b_c = 0; %Ydre radius af kabel
 ElasticModolusCabel = 11.7*10^10; % Elastic Modolus af kobber
 rho_c = 8960; %kg/m3 af kobber
 
%Calculated parameters
lambda = v_0Oil/10^6; % Bølgelængde v/f ved 10 MHz
v_0 = sqrt(c33D/rho_P); % wave speed of compressional waves in the piezoelectric plate
v_0m = sqrt(ElasticModolusMellem/rho_m); % Wave speed of compressional waves in the mellemlag.. nice danlish
v_0b = sqrt(ElasticModolusBacking/rho_b); % Wave speed of compressional waves in the backing material
v_0f = sqrt(ElasticModolusFront/rho_f); % Wave speed of compressional waves in the front layer
v_0a = sqrt(ElasticModolusAdhesive/rho_a); % Wave speed of compressional waves in the adhesive
v_0c = sqrt(ElasticModolusCabel/rho_c); % Wave speed of compressional waves in the adhesive
n=h_33*C_0; % ElectricMatrix parameter.
Z0a = rho_P*v_0*Area; % plane wave acoustic impedance of the piezoelectric plate.
Zba = rho_b*v_0b*Area; % <- vides ikke endnu, men det er: Corresponding acoustic impedance of the backing Z_b^A Vi gætter foreløbigt på =rho_backing*v_0,backing*S Problemet er at der står den er afhængig af frekvens. Plus den er cone-shaped.
Zma = rho_m*v_0m*Area;
Zaa = rho_a*v_0a*Area;
Z0a_f = rho_f*v_0f*Area; % Acoustic impedance front layer ;
ZrAa =v_0Oil*rho_oil*Area; % Acoustic impedance of radiating medium. In this case the oil
kappa_l=1/(rho_oil*v_0Oil^2); %Gorkov lort
kappa_p = 1/(rho_p*v_0p^2); % mere gorkov lort
R_reflect = (rho_f*v_0f-rho_oil*v_0Oil) / (rho_oil*v_0Oil + rho_f*v_0f); % Reflection coefficient between transducer face and oil.
f_1 = 1 - kappa_p/kappa_l; % En konstant der skal bruges i gorkov
f_2 = (2*(rho_p-rho_oil)) / (2*rho_p + rho_oil); % Også en konstant der skal bruges i gorkov
Contrast_factor = (5*rho_p-2*rho_oil) / (2*rho_p+rho_oil) - kappa_p/kappa_l ; 

 %% Parameter der ikke bliver brugt
 %Elastic compliance coefficients [10^-12m^2/N]
%S11T = 11.8*10^(-12);
%        %Frequency coefficients [Hz*m]
%     N_p = 2270;
%     N_1 = 1640;
%     N_3 = 2010;
%     N_t = 2110;
%         % Coupling factors
%     k_p = 0.56;
%     k_t = 0.46;
%     k_31 = 0.32;
%     k_33 = 0.66;
%     k_15 = 0.63;
%         %Piezoelectric voltage Coefficients [10^-3Vm/N]
%     g_31 = -11.2*10^(-3);
%     g_33 = 25*10^(-3);
%      % Piezoelectric Charge Coefficients [10^-12C/N)
%     d_31 = -120*10^(-12);
%     d_15 = 475*10^(-12);
%         %Elastic compliance coefficients [10^-12m^2/N]
%     Q_m = 2000; % Mechanical quality factor
%DielectricLossFactor = 3; % -||- i [10^(-3)]
%V_particle=(4/3)*pi*(r_particle)^3; %Volumen af partikel, bruges i gorkov,
%burde også regnes i gorkov, så den er kommenteret ud her
