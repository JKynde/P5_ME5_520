%% Initialisering
clc
clear
close all
% Det her scribt tjekker for input parameter sensitivitet 
param=StructCreator();
Difftype = "forward"; % Vælg diff type. central, forward eller backward
fields = ["l_aa";"d_33";"l_a";"C_0";"ElasticModolusAdhesive";"l_c";"l_m";"d_p"];
f = 10^6;
V_in = 150;
param.mode=4;
param.Tlmode = 1;
OGparam=param;
BaseCase = abs(Matricer2(f,V_in,param));
scalingfactor = 0.2;
upscale = 1+scalingfactor;
downscale = 1-scalingfactor;
format long
%% Parametervariation

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






