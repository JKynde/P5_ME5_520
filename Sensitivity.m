% Det her scribt tjekker for input parameter sensitivitet
clc
clear
%close all
 %% Running options
param=StructCreator();
Difftype = "forward"; % Vælg diff type. central, forward eller backward
Recalc = 1; % Hvis 0 så laver der bare om på parametrene uden at regne alle udregnede konstanter. 
% Hvis 1, så regner den det hele igennem igen.
%fields = ["l_aa";"d_33";"l_a";"C_0";"ElasticModolusAdhesive";"l_c";"l_m";"d_p";"r_transducer"];
    % Hvis man bare vil have dem alle sammen
    fields=fieldnames(param);
    fields=string(fields);
    fields(1:2)=[]; %Fjerner param.mode og param.Tlmode fra variationen.
scalingfactor = 0.2;
f = 10^6;
V_in = 150;
param.mode=4;
param.Tlmode = 1;

%% Initialisering
OGparam=param;
BaseCase = abs(Matricer2(f,V_in,param));
upscale = 1+scalingfactor;
downscale = 1-scalingfactor;
Resultsup=zeros(1,length(fields));
Resultsdown=zeros(1,length(fields));
Central=zeros(1,length(fields));
Forward=zeros(1,length(fields));
Backward=zeros(1,length(fields));
stepsize=zeros(1,length(fields));
format long
%% Parametervariation
if Recalc == 0
    for n=1:length(fields)
        fprintf("%s is being varied. Its value is originally %e\n",fields(n),param.(fields(n)))
        param.(fields(n))=OGparam.(fields(n))*upscale;
        fprintf("Upping %s with a factor %f. Its value is now %e\n",fields(n),upscale,param.(fields(n)))
        fprintf("Evaluating resulting force output\n")
        Resultsup(n) = abs(Matricer2(f,V_in,param));
        param.(fields(n))=OGparam.(fields(n))*downscale;
        fprintf("Downing %s with a factor %f Its value is now %e \n",fields(n),downscale,param.(fields(n)))
        fprintf("Evaluating resulting force output\n")
        Resultsdown(n)=abs(Matricer2(f,V_in,param));
        stepsize(n)=abs( OGparam.(fields(n))-param.(fields(n)) );
        param=OGparam;
    end
elseif Recalc==1
     for n=1:length(fields)
        fprintf("%s is being varied. Its value is originally %e\n",fields(n),param.(fields(n)))
        param.(fields(n))=OGparam.(fields(n))*upscale;
        fprintf("Upping %s with a factor %f. Its value is now %e\n",fields(n),upscale,param.(fields(n)))
        fprintf("Recalculating parameters")
        param=StructRecalculator(param);
        fprintf("Evaluating resulting force output\n")
        Resultsup(n) = abs(Matricer2(f,V_in,param));
        param.(fields(n))=OGparam.(fields(n))*downscale;
        fprintf("Downing %s with a factor %f Its value is now %e \n",fields(n),downscale,param.(fields(n)))
        fprintf("Recalculating parameters")
        param=StructRecalculator(param);
        fprintf("Evaluating resulting force output\n")
        Resultsdown(n)=abs(Matricer2(f,V_in,param));
        stepsize(n)=abs( OGparam.(fields(n))-param.(fields(n)) );
        param=OGparam;
    end

end
%% Plotting og post processing
for n=1:length(fields)
    Central(n) = abs( (Resultsup(n)-Resultsdown(n))/(2*stepsize(n)) * OGparam.(fields(n))/BaseCase ) ;
    Forward(n) = abs( (Resultsup(n)-BaseCase)/(stepsize(n)) * OGparam.(fields(n))/BaseCase ) ;
    Backward(n) = abs( (BaseCase-Resultsdown(n))/(stepsize(n)) * OGparam.(fields(n))/BaseCase ) ;
end
if Difftype=="central"
    plot=Central;
elseif Difftype=="forward"
    plot=Forward;
elseif Difftype=="backward"
    plot=Backward;
end
%plot=plot/max(plot);
figure;
X = categorical(fields);
X = reordercats(X,fields);
bar(X,plot)
figure;
Bruh(1:length(fields))=BaseCase;
Y = [Resultsup;Bruh;Resultsdown];
Y=Y.';
bar(X,Y)






