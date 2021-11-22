clc %hi hi
clear
fprintf('Initializing \n')
%Parameters;
param=StructCreator();
%% Parameters og running options
plotte = 4; % plotting option 1 er plot mellem F_AC og P-felt til t = 0 - 2 er plot af F_AC_ex og F_AC - 3 er normaliseret plot til at se faserne
generate_F_AC_ex_and_syms = 0; % Vælg om scriptet laver en function til kraften eller om den skal hoppe det over.
Tlmode=1; % Front lag eller ejjjj. Vi elsker voooores børn,
f=10^6; % Frekvens
V_in=150; % Spændingsfaser
omega=2*pi*f; % Vinkelhastighed
t_steps=50; % antal inddelinger i tid, skal bruges til num.int.
param.r_particle = 100*10^-6; % Particle radius should be r<<lambda
z_stepsize=1.2500e-05; %Stepsize til z.
%% Initialisering
V_particle=(4/3)*pi*(param.r_particle)^3; %Volumen af partikel
[F,v_t,~]=Matricer2(f,V_in,param); %Henter F og v fra matricer
z=[param.r_transducer*2:z_stepsize:param.r_transducer*2+2*param.lambda]; %initaliserer afstandsvektoren fra 2 gange transducer radius, til 2 bølgelængder væk
t_stepsize=10^-6/(2*t_steps); %Størrelsen af steppet i tiden til simpsons rule i integral boy.
t=[0:t_stepsize:10^(-6)-t_stepsize]; % initialize time vector t_start:t_step:t_end t løber fra 0 til T (10^-6 s) minus en enkelt stepsize så det passer med at simpson får et lige antal punkter.
U_AC_V=zeros(1,length(z)); %Initialize primary radiation potential per volume vector 
F_AC_V=zeros(1,length(z)); %Initialize force per volume on a particle vector
F_AC=zeros(1,length(z)); %Initialize Force on a particle vector
WavesumresP=zeros(length(t),length(z)); %Laver et array med alle værdier for trykfeltet i både tid og afstand
Wavesumresv=zeros(length(t),length(z)); %Laver et array med alle værdier for hastighedsfeltet i både tid og afstand
WavesumresP_squared=zeros(length(t),length(z)); % Til p^2
Wavesumresv_squared=zeros(length(t),length(z)); % Til v^2
lengthz=length(z);
lengtht=length(t);
%% Long ass nested forloop
Progress=0;
fprintf('Filling wavesumres and squaring \n')
for n=1:lengthz
    for m=1:lengtht
    [~, ~, WavesumresP(m,n),~, ~, ~,Wavesumresv(m,n)]=Pressure(z(n),omega,t(m),F,v_t,param);  
    WavesumresP_squared(m,n)=(WavesumresP(m,n))^2;
    Wavesumresv_squared(m,n)=(Wavesumresv(m,n))^2;
    end
    Progress = (n/(lengthz)*100);
    fprintf('Filling and squaring progress: %f \n',Progress)
end

%% Integration og differentiation
P_avg=zeros(1,lengthz);
v_avg=zeros(1,lengthz);
fprintf('Integrating to find <p^2> and <v^2> \n')

% Integration gøres med simpsons rule
for n=1:lengthz
   P_avg(n) = (1/(1/f))*(integralboy(WavesumresP_squared(:,n),t_stepsize)); % trykket som funktion af tid og afstand i anden, integreret over tiden, divideret med perioden <p^2>
   v_avg(n) = (1/(1/f))*(integralboy(Wavesumresv_squared(:,n),t_stepsize)); % hastigheden som funktion af tid og afstand i anden, integreret over tiden, divideret med perioden <v^2>
end

fprintf('Calculating U_AC_V \n')
for n=1:lengthz
   U_AC_V(n) = (param.f_1/(2*param.rho_oil*(param.v_0Oil)^2)*P_avg(n)-param.f_2*(3/4)*param.rho_oil*v_avg(n)); % Her udregnes gorkovs potential pr volumen af partikel
end

fprintf('Differentiating \n')
F_AC_V=-(differentialboy(U_AC_V,z_stepsize));
F_AC = V_particle*F_AC_V;

%Udregning af hastighed i kvasistatisk betragtning. Det vil sige, at det
%korte transiente forløb hvor partiklen kommer op til steady state farten
%ses der bort fra
%Kommer fra F-Bxdot = m*a og så kigger vi i det scenerie hvor a=0. Dvs
%xdot=F/B og så antager vi at stokes lov gælder.
Kvasihastighed = F_AC/(6*pi*param.mu_oil*param.r_particle); %B fra stokes lov. 

%% Plotting omkring F_AC
if plotte == 1
    fprintf('Plotting \n')
    hold on 
    xlabel('Distance from transducer head')
    yyaxis left
    ylabel('Force on particle with given volume')
    %plot(z,U_AC_V)
    plot(z,Kvasihastighed)
    yyaxis right 
    ylabel('Pressurefield at time = 0')
    %plot(z,P_avg)
    plot(z,WavesumresP(1,:)) 
    %plot (z,U_AC_V)
end

%% F_AC udtryksgenerator - Sinusbølge
if generate_F_AC_ex_and_syms == 1
    % % F_AC function generator. F_AC antages at være en sinus funktion med
    % given amplitude, bølgelængde og fasedrej, som findes ud fra det numeriske plot.
    fprintf('Creating F_AC_function from numeric plot \n')
    zeromap = zeros(1,3); % Laver vektor til n-værdier før fortegnskifte.
    F_AC_max = max(F_AC); % Amplituden af F_AC findes.
    % For loop to check for 3 sign flips (equal to three zeroes) and then
    % lambda can be determined.
    i = 1;
    for n=1:length(F_AC) % Det her for loop går gennem F_AC vektoren og leder efter 
        % fortegnskift. Når den har fundet et fortegnskift på næste indgang,
        % gennem den tilsvarende n-værdi i zeromap og i tæller op. Når den har
        % fundet 3 nuller stopper den med at lede.

        if n ~= length(F_AC) && i <= 3
             if sign(F_AC(n)) ~= sign(F_AC(n+1))
                zeromap(i) = n;
                i=i+1;
             else
             end
        end
    end

    F_AC_lambda = z(zeromap(3))-z(zeromap(1)); % Bølgelængden udregnes, som forskellen mellem 
    % det første og sidste nul, da man har 3 nuller in en enkelt bølgelængde.

    F_AC_ex = F_AC_max*sin((2*pi/F_AC_lambda)*z); % Nu bygges der et udtryk for 
    % F_AC (ex for expression), men den mangler fasedrejet. Bemærk at z
    % vektoren indgør for at F_AC_ex udregnes for de samme z-værdier som F_AC.
    % Derved passer n-værdierne og afstande og funktionsværdier kan
    % sammenlignes.

    zeromap2 = zeros(1,3);
    i = 1;
    for n=1:length(F_AC_ex) % Samme sang og dans gøres nu for F_AC_ex. Altså nullerne findes.

        if n ~= length(F_AC_ex) && i <= 3
             if sign(F_AC_ex(n)) ~= sign(F_AC_ex(n+1))
                zeromap2(i) = n;
                i=i+1;
             else
             end
        end
    end
    deltaZ=z(zeromap(1))-z(zeromap2(1)); % Nu kan vi finde afstanden mellem nullerne og 
    % derved udrenge fasedrejet.
    F_AC_ex = F_AC_max*sin((2*pi/F_AC_lambda)*z-2*pi/F_AC_lambda * deltaZ); % F_AC_ex opdateres
    % med fasedrejet som 2pi/lambda * deltaZ.
    syms x
    F_AC_syms = F_AC_max*sin((2*pi/F_AC_lambda)*x-2*pi/F_AC_lambda * deltaZ); % Symbolic
    % udtryk for F_AC er ovenstående. 
    if plotte ==2 
        plot(z,F_AC_ex,z,F_AC); % Plot til at sammenligne F_AC_ex og F_AC
        %vektorene
    end
end
%% Normaliseret plottefis.
if plotte == 3 || plotte == 4
    fprintf('Normalizing and plotting \n')
    P_avg_normalized=P_avg/max(P_avg);
    v_avg_normalized=v_avg/max(v_avg);
    U_AC_V_normalized=U_AC_V/(min(U_AC_V));  %min fordi U_AC_V er negativ
    F_AC_V_normalized=(F_AC_V) / (max(F_AC_V)); 
    F_AC_normalized=(F_AC) / (max(F_AC));
    
    WavesumresP_normalized=zeros(length(t),length(z));
    Wavesumresv_normalized=zeros(length(t),length(z)); 
    WavesumresP_squared_normalized=zeros(length(t),length(z)); 
    Wavesumresv_squared_normalized=zeros(length(t),length(z)); 
    
    WavesumresP_max=0;
    Wavesumresv_max=0;
    WavesumresP_squared_max=0;
    Wavesumresv_squared_max=0;
    
    for n=1:lengthz
        for m=1:lengtht
           if WavesumresP(m,n) > WavesumresP_max
               WavesumresP_max = WavesumresP(m,n);
           end
           if Wavesumresv(m,n) > Wavesumresv_max
               Wavesumresv_max = Wavesumresv(m,n);
           end
           if WavesumresP_squared(m,n) > WavesumresP_squared_max
               WavesumresP_squared_max = WavesumresP_squared(m,n);
           end
           if Wavesumresv_squared(m,n) > Wavesumresv_squared_max
               Wavesumresv_squared_max = Wavesumresv_squared(m,n);
           end
        end
    end

    for n=1:lengthz
        for m=1:lengtht
            WavesumresP_normalized(m,n)=WavesumresP(m,n)/WavesumresP_max(1,1); 
            Wavesumresv_normalized(m,n)=Wavesumresv(m,n)/Wavesumresv_max; 
            WavesumresP_squared_normalized(m,n)=WavesumresP_squared(m,n)/WavesumresP_squared_max; 
            Wavesumresv_squared_normalized(m,n)=Wavesumresv_squared(m,n)/Wavesumresv_squared_max;
        end
    end
    if plotte == 3
    for n=1:lengtht
        hold on
        axis([2*param.r_transducer  z(length(z)) -1 1]);
        plot(z,WavesumresP_normalized(n,:),'-')
        plot(z,Wavesumresv_normalized(n,:),'--')
        plot(z,P_avg_normalized,'-')
        plot(z,v_avg_normalized,'--')
        plot(z,U_AC_V_normalized,':')  
        plot(z,F_AC_V_normalized,'-.')
        %plot(z,WavesumresP_squared_normalized(n,:),':')
        %plot(z,Wavesumresv_squared_normalized(n,:),'-.')
        if n ~= lengtht
            pause (0.1)
            clf
        else
        end
        end
    end
 end

%% Forsøg på konturplots
if plotte == 4
fprintf('Making pancakes... \n')
%initialisering af plots
[X, Y]=meshgrid(param.r_transducer*2:z_stepsize:param.r_transducer*2+2*param.lambda, 0:1); %her laves meshgrid. Basically den overflade funktionen skal plottes på. Alle inddellinger i længderetningen og to i tidsretningen, da den skal være mindst 2X2
Q = tiledlayout(5,1); %hvor mange grafer man vil have
Q.TileSpacing = 'compact';
title(Q,'Phase Comparison Of Different Quantities', 'fontweight', 'Bold', 'fontsize', 10)
xlabel(Q,'Distance from transducer head', 'fontweight', 'bold', 'fontsize', 10) %samlet label for x-aksen

%kontur af trykfeltet til t=0
ax = nexttile;
pressurefieldkontur=zeros(239,2);
pressurefieldkontur(:,1)=WavesumresP_normalized(1,:);
pressurefieldkontur(:,2)=WavesumresP_normalized(1,:);
contourf(X, Y, pressurefieldkontur', 10); colormap jet;
title('Pressurefield at t=0')
%axis off

%kontur af P_avg
ax1 = nexttile;
peepee=zeros(239,2);
peepee(:,1)=P_avg_normalized(1,:);
peepee(:,2)=P_avg_normalized(1,:);
contourf(X, Y, peepee', 10); colormap jet;
title('Mean square pressure amplitude(normalized)')
%axis off

%kontur af v_avg
ax2 = nexttile;
weewee=zeros(239,2);
weewee(:,1)=v_avg_normalized(1,:);
weewee(:,2)=v_avg_normalized(1,:);
contourf(X, Y, weewee', 10); colormap jet;
title('Mean square velocity amplitude(normalized)')
%axis off

%kontur af U_AC
ax3 = nexttile;
u=zeros(239,2);
u(:,1) = U_AC_V_normalized(1,:);
u(:,2) = U_AC_V_normalized(1,:);
contourf(X, Y, u', 10); colormap jet;
title('Primary radiation potential(normalized)')
%axis off

%kontur af F_particle
ax4 = nexttile;
dongwang=zeros(239,2);
dongwang(:,1) = F_AC_normalized(1,:);
dongwang(:,2) = F_AC_normalized(1,:);
contourf(X, Y, dongwang', 10); colormap jet;
title('Primary radiation force(normalized)');
%axis off
%Kosmetisk grafbehandling
cb = colorbar; %laver en colorbar
cb.Layout.Tile = 'east'; % colorbaren skal være til højre
linkaxes([ax, ax1, ax2, ax3, ax4], 'xy'); %aksene sættes sammen
xticklabels(ax,{}) %her fjernes x og y-aksen for alle grafer, bortset fra x aksen af den nederste
yticklabels(ax,{}) %det er nok den mest komplicerede måde at gøre det på men fuck det
xticklabels(ax1,{})
yticklabels(ax1,{})
xticklabels(ax2,{})
yticklabels(ax2,{})
xticklabels(ax3,{})
yticklabels(ax3,{})
yticklabels(ax4,{})
fprintf('Breakfast is ready! \n')
end
%% done
fprintf('Done \n')


% Welcome to da bottom bruh
% lalala

