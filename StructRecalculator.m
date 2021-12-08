function [paramOut] = StructRecalculator(param)


param.Area = (param.r_transducer)^2*pi; % Plate area
param.epsilon33T = param.PermittivityPolarizationInDirection*param.epsilon_0; %epsilon 33T skal bruges senere til at udregnes h_33, som vi skal bruge.
param.h_33=param.d_33/(param.S33E*param.epsilon33T); % Formel for h_33, som vi tror er Piezoelectric stiffness constant er regnet ud fra side 15 i https://link.springer.com/content/pdf/bbm%3A978-94-007-0579-1%2F1.pdf
param.C_0 = param.Area/(1/(param.epsilon33T)*param.d_p);
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
param.Zma = param.rho_m*param.v_0m*param.Area;
param.Zaa = param.rho_a*param.v_0a*param.Area;
param.Z0a_f = param.rho_f*param.v_0f*param.Area; % Acoustic impedance front layer ;
param.ZrAa =param.v_0Oil*param.rho_oil*param.Area; % Acoustic impedance of radiating medium. In this case the oil

param.kappa_l=1/(param.rho_oil*param.v_0Oil^2); %Gorkov lort
param.kappa_p = 1/(param.rho_p*param.v_0p^2); % mere gorkov lort

param.R_reflect = (param.rho_f*param.v_0f-param.rho_oil*param.v_0Oil) / (param.rho_oil*param.v_0Oil + param.rho_f*param.v_0f); % Reflection coefficient between transducer face and oil.

param.f_1 = 1 - param.kappa_p/param.kappa_l; % En konstant der skal bruges i gorkov
param.f_2 = (2*(param.rho_p-param.rho_oil)) / (2*param.rho_p + param.rho_oil); % Også en konstant der skal bruges i gorkov

param.Contrast_factor = (5*param.rho_p-2*param.rho_oil) / (2*param.rho_p+param.rho_oil) - param.kappa_p/param.kappa_l; 

param.mu_oil = param.nu_oil*param.rho_oil;

paramOut=param;
end

