clc
clear
OGparam=StructCreator(); % Originale parametre indlæses
%% Running options
V_in=150;
f=1.14*10^6; % Creates a linear space of f values between
OGparam.mode=3; % Indstil running mode for matricer 2
OGparam.Tlmode=1; % Front lag eller eeeejjj.
halvesmax=100; % Maksimale antal halveringer af stepsize
itermax=2;
scaling = 0.2; % Initial scaling factor to determine initial stepsize for each parameter
Bounds = 0;
Boundscale = 2; % Scaleringsfaktor hvormed bounds defineres
fields = ["l_aa";"l_a";"l_m";'d_p']; % Parametre som varieres
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
BaseCaseResult=abs(Matricer2(f,V_in,CurrentGuess)); %Samme for resultaterne
BestGuessResult=BaseCaseResult;
CurrentGuessResult=BaseCaseResult;
halves=0; % Antal halveringer af stepsize
iter=0;
%% Big ol' while loop
while stop==0
    iter=iter + 1;
    for n=1:lengthfields % For loopet gennemløber alle valgte fields.
        CurrentGuess.(fields(n))=CurrentGuess.(fields(n))+stepsize(n); %Feltet steppes op med stepsize
        if Bounds==1 % Hvis bound er til så limit 
            if CurrentGuess.(fields(n))>Upperbound(n)
                CurrentGuess.(fields(n))=Upperbound(n);
            end
        end
        CurrentGuess=StructRecalculator(CurrentGuess); %Recalc params
        CurrentGuessResult=abs(Matricer2(f,V_in,CurrentGuess)); % Evaluer.
            if CurrentGuessResult>BestGuessResult
                BestGuess=CurrentGuess;
                BestGuessResult=CurrentGuessResult;
            end
        CurrentGuess=BaseCase;
        CurrentGuessResult=BaseCaseResult;
    end
    CurrentGuess=BaseCase;
    CurrentGuessResult=BaseCaseResult;
    for n=1:lengthfields
        CurrentGuess.(fields(n))=CurrentGuess.(fields(n))-stepsize(n); % Step down
        if Bounds==1
            if CurrentGuess.(fields(n))<Lowerbound(n)
                CurrentGuess.(fields(n))=Lowerbound(n);
            end
        end
        CurrentGuess=StructRecalculator(CurrentGuess); % Recalc
        CurrentGuessResult=abs(Matricer2(f,V_in,CurrentGuess)); %Evaluer
        if CurrentGuessResult>BestGuessResult
            BestGuess=CurrentGuess;
            BestGuessResult=CurrentGuessResult;
        end
        CurrentGuess=BaseCase;
        CurrentGuessResult=BaseCaseResult;
    end

    if BaseCaseResult==BestGuessResult
        stepsize=stepsize./2;
        halves=halves+1;
    end
    if halves>=halvesmax
        stop=1;
    end
    if iter>=itermax
stop=1;
    end
     BaseCase=BestGuess;
     BaseCaseResult=BestGuessResult;
end

for n=1:lengthfields
    if BestGuess.(fields(n))~=OGparam.(fields(n))
        fprintf('A changed Value was %s, which originally was %e and now is %e\n',fields(n),OGparam.(fields(n)),BestGuess.(fields(n)))
    end

end
OGresult=abs(Matricer2(f,V_in,OGparam));
fprintf('The original result was %f and the new result is %f\n',OGresult,BestGuessResult)

