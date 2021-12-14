clc
clear
OGparam=StructCreator(); % Originale parametre indlæses
% OGparam.r_tranducer=0.00485/2;
% OGparam.d_p=0.00205;
% OGparam.ElasticModolusAdhesive=5.9*10^9;
% OGparam.l_aa=0.025*10^-3;
load("ExImpedance.mat");
f_ex=Impedance(:,1);
Z_ex=Impedance(:,2)./Impedance(:,3);
%% Running options
V_in=150;
ZfitorF=1;
if ZfitorF==1
    f=f_ex;
else
    f=1.124890e+06; % Creates a linear space of f values between
end
OGparam.mode=3; % Indstil running mode for matricer 2
OGparam.Tlmode=1; % Front lag eller eeeejjj.
halvesmax=1000; % Maksimale antal halveringer af stepsize
itermax=100; % Maksimale antal iterationer
scaling = 0.3; % Initial scaling factor to determine initial stepsize for each parameter
Bounds = 1; % Bounds?
dontgoneg = 1; % Prevent negative parameters

Boundscale = 3; % Scaleringsfaktor hvormed bounds defineres
fields = ["l_aa";"l_a";"l_m";"ElasticModolusBacking";"rho_b";"rho_f";"ElasticModolusFront";"rho_a";"ElasticModolusAdhesive";"rho_m";"ElasticModolusMellem";"C_0";"d_p";"rho_P";"d_33";"epsilon33T";"c33D";"h_33";"S33E"]; % Parametre som varieres
% fields=fieldnames(OGparam);
% fields=string(fields);
% fields(1:2)=[]; %Fjerner param.mode og param.Tlmode fra variationen.
%% Initialization
lengthfields = length(fields); %
if Bounds==1
    Upperbound = zeros(1,length(fields));
    Lowerbound = zeros(1,length(fields));
    for n=1:lengthfields
        Upperbound(n)=OGparam.(fields(n))*Boundscale; % Øvre og lower bounds laves (Bruge pt ikke)
        Lowerbound(n)=OGparam.(fields(n))./Boundscale;
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
        if dontgoneg==1
            if CurrentGuess.(fields(n))<0
                CurrentGuess.(fields(n))=0;
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
bruh=1;
for n=1:lengthfields
    if BestGuess.(fields(n))~=OGparam.(fields(n))

        indexchangedvalues(bruh)=n;
        bruh=bruh+1;
        fprintf('A changed Value was %s, which originally was %e and now is %e\n',fields(n),OGparam.(fields(n)),BestGuess.(fields(n)))
    end

end

%% Plotting
%Initialisering af plots
if ZfitorF==1
    F_sim=zeros(1,length(f));
    Z_sim=zeros(1,length(f));
    OGZ_sim=zeros(1,length(f));
    OGresult=Zfit(f,V_in,OGparam,Z_ex);
    for n=1:length(f)
        [~,~,OGZ_sim(n)]=Matricer2(f(n),V_in,OGparam); % henter orignal simuleret impedans
    end
    for n=1:length(f)
        [F_sim(n),~,Z_sim(n)]=Matricer2(f(n),V_in,BestGuess); % nyt simuleret plot
    end
    Z_sim=abs(Z_sim);
    F_sim=abs(F_sim);
    OGZ_sim=abs(OGZ_sim);
    figure;
    plot(f, Z_ex, f, Z_sim, f, OGZ_sim);
    title('Result of impedance amplitude curve fitting', 'FontSize', 20)
    ylabel('Impedance [\Omega]', 'FontSize', 20)
    xlabel('Frequency [Hz]', 'FontSize', 20)
    legend({'Measured impedance', 'New simulated impedance', 'Original simulated impedance'}, 'FontSize', 14)
    xlim([10^5 1.2*10^6]);
    xticks([10^5:10^5:1.2*10^6]);


    f_b = logspace(5,6.5,5000);
    for n=1:length(f_b)
        [F_sim(n),~,Z_sim(n)]=Matricer2(f_b(n),V_in,BestGuess);
    end
    F_db = 20*log10(abs(F_sim));
    AngelFoutt = unwrap(rad2deg(angle(F_sim)));
    figure;

    t=tiledlayout(2,1); % Her plottes outputtet af matricer som bodeplot
    xlabel(t,'log spaced f values')

    nexttile
    semilogx(f_b,F_db) % skriv enten F_DB, V_DB eller Z_DB
    grid on
    title('Magnitude plot of F(f)')
    ylabel('Magnitude of F(dB)')

    nexttile
    semilogx(f_b,AngelFoutt)
    grid on
    title('Angle plot of F(f)')
    ylabel('Angle of F(^\circ)')
else
    OGresult=abs(Matricer2(f,V_in,OGparam));
end
fprintf('The original result was %f and the new result is %f\n',OGresult,BestGuessResult)

for n=1:length(indexchangedvalues)
    changedfields(n)=fields(indexchangedvalues(n));
end

figure
X=categorical(changedfields);
X=reordercats(X,changedfields);
for n=1:length(changedfields)
    Y(1,n)=OGparam.(changedfields(n)) / OGparam.(changedfields(n));
    Y(2,n)=BestGuess.(changedfields(n))/OGparam.(changedfields(n));
end
bar(X,Y');
%set(gca,'XtickLabel',changedfields)
title('Original parameters compared to changed parameters', 'FontSize', 20)
ylabel('Normalized parameter values', 'FontSize', 18)
legend('Original parameter', 'changed parameter')





