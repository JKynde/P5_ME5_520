%Definitions
j = sqrt(-1);

% Priezoelektriske Materiale konstanter (PIC181).
rho_P = 7800; % Densitet i [kg/m^3]
epsilon_0 = 8.854*10^(-12); % permittivity of free space [F/m]
PermittivityPolarizationInDirection = 1200; % Relative permittivity in the polarization direction = epsilon_33^T/epsilon_0
epsilon33T = PermittivityPolarizationInDirection*epsilon_0; %epsilon 33T skal bruges senere til at udregnes h_33, som vi skal bruge.
PermittivityNormalDirection = 1500; % Relative permittivity in the norm direction epsilon_11^T/epsilon_0
DielectricLossFactor = 3; % -||- i [10^(-3)]

% Coupling factors
k_p = 0.56;
k_t = 0.46;
k_31 = 0.32;
k_33 = 0.66;
k_15 = 0.63;

% Piezoelectric Charge Coefficients [10^-12C/N)
d_31 = -120*10^(-12);
d_33 = 265*10^(-12);
d_15 = 475*10^(-12);

%Piezoelectric voltage Coefficients [10^-3Vm/N]
g_31 = -11.2*10^(-3);
g_33 = 25*10^(-3);

    % Acousto-mechenical properties
%Frequency coefficients [Hz*m]
N_p = 2270;
N_1 = 1640;
N_3 = 2010;
N_t = 2110;

%Elastic compliance coefficients [10^-12m^2/N]
S11T = 11.8*10^(-12);
S33E = 14.2*10^(-12);
c33D = 16.6*10^10; % Elastic stiffness coefficient [10^10N/m^2]
Q_m = 2000; % Mechanical quality factor
C_0 = 0.104*10^(-9); % Kapacitans af pladen i nanofarad. Det er to af dem så det noget fuck endnu
h_33=d_33/(S33E*epsilon33T); % Formel for h_33, som vi tror er Piezoelectric stiffness constant er regnet ud fra side 15 i https://link.springer.com/content/pdf/bbm%3A978-94-007-0579-1%2F1.pdf
% Der er en antagelse at h_33 i Ultrasonic Nondestructive Evaluation
% systems er den samme som i den artiklen hvor formlen for h_33 er fundet.

%Geometri af en enkelt piezoplade(vi har to. Ved ikke hvordan vi lige
%modellerer det endnu)
d_p = 0.002; %Plate thickness [m]
r_transducer = 0.005/2;
Area = (r_transducer)^2*pi; % Plate area

%Backing
ElasticModolusBacking = 2.1*10^11; %Backing materiale konstanter.[N/m^2]
rho_b = 7800;%Backing densitet [kg/M^3]m

%Front layer parameters. Altså housingen ude foran.
l_a = 3.2*10^(-3); %layer thickness [m]
rho_f = 2700 ; %Densitet af front layer. [kg/m^3]
ElasticModolusFront = 6.9*10^10; %Young's modolus for front end materialet'

%Partikle parameters: sat til luft for sjov for nu
rho_p = 1.225; %kg/m^3
v_0p = 340.19744;% m/s eller 761 i retard units  

% Oil parameters
rho_oil = 857; % Densitet af olien [kg/m^3]
v_0Oil = 1493; % Lydens hastighed i olien.

%Calculated parameters
lambda = v_0Oil/10^6; % Bølgelængde v/f ved 10 MHz
v_0 = sqrt(c33D/rho_P); % wave speed of compressional waves in the piezoelectric plate
v_0b = sqrt(ElasticModolusBacking/rho_b); % Wave speed of compressional waves in the backing material
v_0f = sqrt(ElasticModolusFront/rho_f); % Wave speed of compressional waves in the front layer
n=h_33*C_0; % ElectricMatrix parameter.
Z0a = rho_P*v_0*Area; % plane wave acoustic impedance of the piezoelectric plate.
Zba = rho_b*v_0b*Area; % <- vides ikke endnu, men det er: Corresponding acoustic impedance of the backing Z_b^A Vi gætter foreløbigt på =rho_backing*v_0,backing*S Problemet er at der står den er afhængig af frekvens. Plus den er cone-shaped.
Z0a_f = rho_f*v_0f*Area; % Acoustic impedance front layer ;
ZrAa =v_0Oil*rho_oil*Area; % Acoustic impedance of radiating medium. In this case the oil
kappa_l=1/(rho_oil*v_0Oil^2); %Gorkov lort
kappa_p = 1/(rho_p*v_0p^2);