% Det her scribt tjekker for input parameter sensitivitet
clc
clear
%close all
 %% Running options
param=StructCreator();
Difftype = "central"; % Vælg diff type. central, forward eller backward
Recalc = 1; % Hvis 0 så laver der bare om på parametrene uden at regne alle udregnede konstanter. 
% Hvis 1, så regner den det hele igennem igen.
%fields = ["l_aa";"d_33";"l_a";"C_0";"ElasticModolusAdhesive";"l_c";"l_m";"d_p";"r_transducer"];
    % Hvis man bare vil have dem alle sammen
    fields=fieldnames(param);
    fields=string(fields);
    fields(1:2)=[]; %Fjerner param.mode og param.Tlmode fra variationen.
% fields=["Area"]
scalingfactor = 0.8; % hvor meget hver parameter scaleres med. Skal være 'lille'(definition mangler). Skal være mellem 0 og 1.
f = 1.14*10^6; %resonansfrekvensen
V_in = 150; %input voltage
param.mode=3; % hvor meget der er med i modellen 
param.Tlmode = 1; % med frontlag

%% Initialisering
OGparam=param; % initialisering af parametre, you know the drill
BaseCase = abs(Matricer2(f,V_in,param)); %regner basecase ud fra model
upscale = 1+scalingfactor; %scale up
downscale = 1-scalingfactor; %scale down
Resultsup=zeros(1,length(fields)); %definerer en masse størrelser som nulvektorer
Resultsdown=zeros(1,length(fields));
Central=zeros(1,length(fields));
Forward=zeros(1,length(fields));
Backward=zeros(1,length(fields));
stepsize=zeros(1,length(fields));
format long % lang diller
%% Parametervariation
if Recalc == 0 %uden at rekalkulere. ikke så relevant
    for n=1:length(fields)
        fprintf("%s is being varied. Its value is originally %e\n",fields(n),param.(fields(n)))
        param.(fields(n))=OGparam.(fields(n))*upscale;
        fprintf("Upping %s with a factor %f. Its value is now %e\n",fields(n),upscale,param.(fields(n)))
        fprintf("Evaluating resulting force output\n")
        Resultsup(n) = abs(Matricer2(f,V_in,param));
        param.(fields(n))=OGparam.(fields(n))*downscale;
        fprintf("Downing %s with a factor %f Its value is now %e \n",fields(n),downscale,param.(fields(n)))
        stepsize(n)=abs( OGparam.(fields(n))-param.(fields(n)) );
        fprintf("Evaluating resulting force output\n")
        Resultsdown(n)=abs(Matricer2(f,V_in,param));
        param=OGparam;
    end
elseif Recalc==1 % med rekalkulering
     for n=1:length(fields)
        fprintf("%s is being varied. Its value is originally %e\n",fields(n),param.(fields(n)))
        param.(fields(n))=OGparam.(fields(n))*upscale; %ganger den originale parameter med opskaleringsfaktoren
        fprintf("Upping %s with a factor %f. Its value is now %e\n",fields(n),upscale,param.(fields(n)))
        fprintf("Recalculating parameters")
        param=StructRecalculator(param); %rekalkulerer udregnede parametre med den nye parameter
        fprintf("Evaluating resulting force output\n")
        Resultsup(n) = abs(Matricer2(f,V_in,param)); %udregner ny F med nye parametre
        param.(fields(n))=OGparam.(fields(n))*downscale; %nedskalerer
        fprintf("Downing %s with a factor %f Its value is now %e \n",fields(n),downscale,param.(fields(n)))
        stepsize(n)=abs( OGparam.(fields(n))-param.(fields(n)) ); %udregner stepsizen
        fprintf("Recalculating parameters")
        param=StructRecalculator(param); %rekalkulerer
        fprintf("Evaluating resulting force output\n")
        Resultsdown(n)=abs(Matricer2(f,V_in,param)); %udregner F med nye parametre
        param=OGparam; %resetter parameter, så er klar til at variere ny parameter
    end

end
%% Plotting og post processing
for n=1:length(fields)
    Central(n) = abs( (Resultsup(n)-Resultsdown(n))/(2*stepsize(n)) * OGparam.(fields(n))/BaseCase ) ;
    Forward(n) = abs( (Resultsup(n)-BaseCase)/(stepsize(n)) * OGparam.(fields(n))/BaseCase ) ;
    Backward(n) = abs( (BaseCase-Resultsdown(n))/(stepsize(n)) * OGparam.(fields(n))/BaseCase ) ;
end
% figure;
% X = categorical(fields);
% X = reordercats(X,fields);
% Bruh(1:length(fields))=BaseCase;
% Y = [Resultsup;Bruh;Resultsdown];
% Y=Y.';
% bar(X,Y)
%% Fjerner nuller
m = 1; % initialisering af m og n
n = 0;
if Difftype == 'central'
for n=1:length(fields)
    if Central(n) == 0
        Delvect(m)=n;
        m=m+1;
    end
end
Central(Delvect) = [];
delfields = fields(Delvect);
n=0;
while 1
    n=n+1;
    for m=1:length(delfields)
        if strcmp(fields(n),delfields(m))
            fields(n) = [];
            n=n-1;
        end
    end
    if n == length(fields)
        break
    end
end
end

if Difftype == 'forward'
for n=1:length(fields)
    if Forward(n) == 0
        Delvect(m)=n;
        m=m+1;
    end
end
Forward(Delvect) = [];
delfields = fields(Delvect);
n = 0;
while 1
    n=n+1;
    for m=1:length(delfields)
        if strcmp(fields(n),delfields(m))
            fields(n) = [];
            n=n-1;
        end
    end
    if n == length(fields)
        break
    end
end
end

if Difftype == 'backward'
for n=1:length(fields)
    if Central(n) == 0
        Delvect(m)=n;
        m=m+1;
    end
end
Backward(Delvect) = [];
delfields = fields(Delvect);
n = 0;
while 1
    n=n+1;
    for m=1:length(delfields)
        if strcmp(fields(n),delfields(m))
            fields(n) = [];
            n=n-1;
        end
    end
    if n == length(fields)
        break
    end
end
end
    
%% plotting
if Difftype=="central"
    plot=Central;
elseif Difftype=="forward"
    plot=Forward;
elseif Difftype=="backward"
    plot=Backward;
end
%plot=plot/max(plot);
figure;
B = sort(plot, "descend");
X = categorical(fields);
X = reordercats(X,fields);
bar(X,B)
ylabel('Sensitivity coefficient \phi_i', 'FontSize', 20)
title('Sensitivity of relevant parameters(central difference, \zeta=0.1)', 'FontSize', 20)
% figure;
% Bruh(1:length(fields))=BaseCase;
% Y = [Resultsup;Bruh;Resultsdown];
% Y=Y.';
% bar(X,Y)






