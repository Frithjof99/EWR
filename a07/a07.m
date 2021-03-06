clear;
clc;
close all;

openfigure(4, 'init');

Stadtteile = struct([geopoint()]);  
for i=1:27
    Stadtteile(i) = struct(gpxread(sprintf('data/data%dgpx.sec',i)));
end

Namen = {'Groetzingen', 'Wolfartsweier', 'Oberreut', ...
    'Gruenwettersbach', 'Rueppur', 'Palmbach', 'Hohenwettersbach', ...
    'Waldstadt','Weststadt','Weiherfeld-Dammerstock','Innenstadt-West',...
    'Nordstadt','Daxlanden','Nordweststadt','Suedstadt','Stupferich',...
    'Muehlburg','Rintheim','Gruenwinkel','Knielingen','Durlach',...
    'Suedweststadt','Oststadt','Beiertheim-Bulach','Innenstadt-Ost',...
    'Neureut','Hagsfeld'};

N = [9138, 3156, 9554, 4082, 10630, 1936, 3023, 12484, 20489, ... 
    6029, 10283,9770,11695, 11755, 20121, 2782, 17149, 5991, 10709, ...
    10137, 30473, 20709, 22808, 6974, 6725, 18877, 7140];



phi = GetPhi(0.95);

N_tot = zeros(27,1);
for i = 1:27
    for m = 1:27
        N_tot(i) = N_tot(i) + phi(m,i)*N(m); 
    end
end


T = 100;
dT = 1;

s = T / dT;

dS = zeros(27,s);
dI = zeros(27,s);
dR = zeros(27,s);

S = zeros(27,s);
I = zeros(27,s);
R = zeros(27,s);

S(:,1) = N;

S(16,1) = S(16,1) - 200;
I(16,1) = 200;

beta = 0.5;
gamma = 0.25;

for i = 1:s
    
    %X = Change(beta, gamma, phi, N_tot, S(:,i), I(:,i), R(:,i));
    
    [dS(:,i),dI(:,i),dR(:,i)] = Change(beta, gamma, phi, N_tot, S(:,i), I(:,i), R(:,i));
    S(:,i+1) = S(:,i) + dT * dS(:,i);
    I(:,i+1) = I(:,i) + dT * dI(:,i);
    R(:,i+1) = R(:,i) + dT * dR(:,i);
    
    
    
end

t = 1:s;

figure(2);
hold on;
plot(t*dT,sum(S(:,t) + I(:,t) + R(:,t)));
plot(t*dT,sum(S(:,t)));
plot(t*dT,sum(I(:,t)));
plot(t*dT,sum(R(:,t)));
hold off;
legend("Gesamtbevölkerung", "Anfällig", "Infiziert", "Erholte", "FontSize", 14);
xlabel("Zeit", "FontSize", 14);
ylabel("Bevölkerung", "FontSize", 14);


print('-f2','bild2-1','-dpng','-r800');


figure(3);
clf;
hold on;
plot(t*dT,S(16,t));
plot(t*dT,I(16,t));
plot(t*dT,R(16,t));
hold off;
legend("Anfällig", "Infiziert", "Erholte", "FontSize", 14);
xlabel("Zeit", "FontSize", 14);
ylabel("Bevölkerung", "FontSize", 14);
title("Stadtteil 16", "FontSize", 14);


print('-f3','bild5-1','-dpng','-r400');

figure(4);
clf;
hold on;
plot(t*dT,S(17,t));
plot(t*dT,I(17,t));
plot(t*dT,R(17,t));
hold off;
legend("Anfällig", "Infiziert", "Erholte", "FontSize", 14);
xlabel("Zeit", "FontSize", 14);
ylabel("Bevölkerung", "FontSize", 14);
title("Stadtteil 17", "FontSize", 14);


print('-f4','bild5-2','-dpng','-r400');

I_max = max(I);

h = [1, 10, 20, 30, 40,50,60,70,80,90,100];

size(h)

for i = 1:size(h,2)
    figure(1);
    clf;
    DrawMapFrame(I(:,h(i)), I_max, Stadtteile);
    title(num2str(h(i)),"FontSize", 20);
    %print('-f1',['bild1-' num2str(h(i))],'-dpng','-r400'); %speichert ein Bild
    %pro Zeitschritt
    
    pause(0.0);
    
end

