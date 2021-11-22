%% Initializing
clc
clear
f=10^6;
V_in=150;
param=StructCreator();
param.l_aa = 0.015*10^(-3); % Initial guess til given parameter. Vi starter med limtykkelsen
Upperbound = 1*10^-2; % Upper bound for given parameter.
Lowerbound = 0; % Lower bound for given parameter.
stepsize = 3*10^-3; % Initial stepsize.
stop = 0;
CurrentGuess = param.l_aa;
CurrentGuessValue = 0;
Stepup = 0;
Stepupvalue = 0;
Stepdown = 0;
StepdownValue = 0;

while stop==0

    CurrentGuessValue = abs(Matricer(f,V_in,param)); % Initial/Current Guess value of F

    Stepup=CurrentGuess + stepsize;


    Stepdown=CurrentGuess + stepsize;





end


























