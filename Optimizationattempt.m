clc
clear
OGparam=StructCreator(); % Originale parametre indlæses
load("ExImpedance.mat");
f_ex=Impedance(:,1);
Z_ex=Impedance(:,2)./Impedance(:,3);
%% Running options
V_in=150;
ZfitorF=1;
if ZfitorF==1
    f=f_ex;
else
    f=1.14*10^6; % Creates a linear space of f values between
end
OGparam.mode=3; % Indstil running mode for matricer 2
OGparam.Tlmode=1; % Front lag eller eeeejjj.
halvesmax=1000; % Maksimale antal halveringer af stepsize
itermax=50;
scaling = 0.1; % Initial scaling factor to determine initial stepsize for each parameter
Bounds = 1;
Boundscale = 2; % Scaleringsfaktor hvormed bounds defineres
%fields = ["l_aa";"l_a";"l_m";'d_p';"rho_P"]; % Parametre som varieres
fields=fieldnames(OGparam);
fields=string(fields);
fields(1:2)=[]; %Fjerner param.mode og param.Tlmode fra variationen.
%% Initialization
lengthfields = length(fields); %
if Bounds==1
    Upperbound = zeros(1,length(fields));
    Lowerbound = zeros(1,length(fields));
    for n=1:lengthfields
        Upperbound(n)=OGparam.(fields(n))*Boundscale; % Øvre og lower bounds laves (Bruge pt ikke)
        Lowerbound(n)=OGparam.(fields(n))/Boundscale;
    end
end
stepsize=zeros(1,lengthfields);
for n=1:lengthfields % Der findes en individuel stepsize til hver parameter
    stepsize(n)=OGparam.(fields(n))*scaling;%-OGparam.(fields(n));
end
stop = 0; % Stop value
BaseCase=OGparam; % Best Case struct initialiseres til OGparam
CurrentGuess=OGparam; % Samme for nuværende gæt
BestGuess=CurrentGuess; % Samme for beste gæt
if ZfitorF==1
    BaseCaseResult=Zfit(f,V_in,CurrentGuess,Z_ex);
elseif ZfitorF == 0
    BaseCaseResult=abs(Matricer2(f,V_in,CurrentGuess)); %Samme for resultaterne
end
BestGuessResult=BaseCaseResult;
CurrentGuessResult=BaseCaseResult;
halves=0; % Antal halveringer af stepsize
iter=0;
%% Big ol' while loop
while stop==0
    iter=iter + 1
    for n=1:lengthfields % For loopet gennemløber alle valgte fields.
        CurrentGuess.(fields(n))=CurrentGuess.(fields(n))+stepsize(n); %Feltet steppes op med stepsize
        if Bounds==1 % Hvis bound=1 er til så er der limits på
            if CurrentGuess.(fields(n))>Upperbound(n) %hvis current gæt er større end bounds, sæt current gæt til upper limit
                CurrentGuess.(fields(n))=Upperbound(n);
            end
        end
        CurrentGuess=StructRecalculator(CurrentGuess); %Recalc params
        if ZfitorF==1
            CurrentGuessResult=Zfit(f,V_in,CurrentGuess,Z_ex); % Evaluer
            if CurrentGuessResult<BestGuessResult
                BestGuess=CurrentGuess;
                BestGuessResult=CurrentGuessResult; %og gem resultat, hvor stor kraften er ved det bedste guess
            end
        else
            CurrentGuessResult=abs(Matricer2(f,V_in,CurrentGuess)); % Evaluer.
            if CurrentGuessResult>BestGuessResult %hvis current gæt er bedre end det hidtil bedste, gem current gæt som det hidtil bedste
                BestGuess=CurrentGuess;
                BestGuessResult=CurrentGuessResult; %og gem resultat, hvor stor kraften er ved det bedste guess
            end
        end
        CurrentGuess=BaseCase;
        CurrentGuessResult=BaseCaseResult;
    end
    CurrentGuess=BaseCase; %reset current gæt til base gæt
    CurrentGuessResult=BaseCaseResult; %reset current gæt resultat til base gæt
    for n=1:lengthfields
        CurrentGuess.(fields(n))=CurrentGuess.(fields(n))-stepsize(n); % Step down
        if Bounds==1 %hvis bounds=1, sæt nedre grænse for current gæt
            if CurrentGuess.(fields(n))<Lowerbound(n)
                CurrentGuess.(fields(n))=Lowerbound(n);
            end
        end
        CurrentGuess=StructRecalculator(CurrentGuess); % Recalc
        if ZfitorF==1
            CurrentGuessResult=Zfit(f,V_in,CurrentGuess,Z_ex); % Evaluer
            if CurrentGuessResult<BestGuessResult
                BestGuess=CurrentGuess;
                BestGuessResult=CurrentGuessResult; %og gem resultat, hvor stor kraften er ved det bedste guess
            end
        else
            CurrentGuessResult=abs(Matricer2(f,V_in,CurrentGuess)); %Evaluer
            if CurrentGuessResult>BestGuessResult %Hvis current gæt result er større end det hidtil største resultat, gem dette som det bedste resultat
                BestGuess=CurrentGuess;
                BestGuessResult=CurrentGuessResult; % og gem selve resultatet så
            end
        end
        CurrentGuess=BaseCase; % current gæt resettes inde i forloopet
        CurrentGuessResult=BaseCaseResult; %current gæt resultat resettes inde i forloopet
    end
    CurrentGuess=BaseCase; %current gæt resettes udenfor forloop
    CurrentGuessResult=BaseCaseResult; %current gæt resultat resettes udenfor forloop

    if BaseCaseResult==BestGuessResult %hvis base er lig med det bedste, dsv ingen nye veje er bedre end der hvor vi er, halver stepsize
        stepsize=stepsize./2;
        halves=halves+1 %tæl halveringer
    end
    if halves>=halvesmax %sæt max for halveringer
        stop=1;
    end
    if iter>=itermax %sæt max for iterationer
        stop=1;
    end
    BaseCase=BestGuess; %gem det bedste resultat fra denne iteration, som den nye basecase
    BaseCaseResult=BestGuessResult;
end

for n=1:lengthfields
    if BestGuess.(fields(n))~=OGparam.(fields(n))
        fprintf('A changed Value was %s, which originally was %e and now is %e\n',fields(n),OGparam.(fields(n)),BestGuess.(fields(n)))
    end

end
if ZfitorF==1
    F_sim=zeros(1,length(f));
    Z_sim=zeros(1,length(f));
    OGresult=Zfit(f,V_in,OGparam,Z_ex);
    for n=1:length(f)
        [F_sim(n),~,Z_sim(n)]=Matricer2(f(n),V_in,BestGuess);
    end
    Z_sim=abs(Z_sim);
    F_sim=abs(F_sim);

else
    OGresult=abs(Matricer2(f,V_in,OGparam));
end
fprintf('The original result was %f and the new result is %f\n',OGresult,BestGuessResult)


