%% Aufgabe5 (Von der Erde zum Mond)
%In dieser Aufgabe wird die Reise von der Erde zum Mond modelliert und
%numerisch berechnet.

%Autor: Lena Hilpp ; Jan Frithjof Fleischhammer
%Version: 11.05.2020

%% Teilaufgabe a
%Reset Workspace, Reset Konsole, alle Fenster schliessen
clear;clc;close all;

%Initialisierung der benoetigten Fenster
openfigure(7,'init'); 

%% Teilaufgabe b
%Konstanten definieren umgerechnung mit Referenzgroessen
mE=1; %Masse Erde in kg
mM=1.23*10^-2; %Masse Mond in kg
RE=1.659*10^-2;  %Erdradius in m
RM=4.521*10^-3;  %Mondradius in m
y=86400^2*(3.844)^(-3)*5.974*6.673*10^(-11); %Gravitationskonstante un m^3/(kg s^2)

%Parameter definieren
v0=86400*1.109/(3.844*10^4);  %Abschussgeschwindigkeit in m/s
theta=0.008;      %Abschusswinkel in rad

%Entfernung Erde
dE=@(x1,x2) sqrt((x1-0).^2+(x2-0).^2);
%Entfernung Mond
dM=@(x1,x2) sqrt((x1-1).^2+(x2-0).^2);
%Gravitationspotential V
pot=@(x1,x2)-y*((mE./dE(x1,x2))+(mM./dM(x1,x2)));
%Gravitationsfeldstaerke G
frc=@(x1,x2)[-y.*(mE./((x1.^2+x2.^2).^(3/2)).*x1)-y.*(mM./(((x1-1).^2+x2.^2).^(3/2)).*(x1-1)),-y.*(mE./((x1.^2+x2.^2).^(3/2)).*x2)-y.*(mM./(((x1-1).^2+x2.^2).^(3/2)).*x2)];

%Plot Gravitationspotential
figure(1);
clf;
Nx=150;
x = linspace(0,1.2,Nx)';
y = linspace(-0.5,0.5,Nx)';
[X,Y]=meshgrid(x,y);
Z=pot(X,Y);
R=dE(X,Y);
Z(R<5*RE)=NaN;
mesh(X,Y,Z);
print('-f1','bild1','-dpng','-r100');
%Plot Gravitationsfeldstaerke
figure(2);
clf;
x = linspace(-0.2,1.2,50)';
y = linspace(-0.5,0.5,40)';
[X,Y]=meshgrid(x,y);
A=frc(X,Y);
R1=dE(X,Y);
A(R1<5*RE)=NaN;
idx=length(A(1,:))/2;
U=A(:,1:idx);
V=A(:,idx+1:end);
quiver(X,Y,U,V);
print('-f2','bild2','-dpng','-r800');
%% Teilaufgabe c
x0=[RE,0]; %Startposition
V0=v0*[cos(theta),sin(theta)]; %Abschussgeschwindigkeit
%Anfangswertproblem
Fode= @(tt,zz) [zz(3);zz(4);frc(zz(1),zz(2))'];
z0=[x0,V0];

%% Teilaufgabe d
%Loesen der Diff.Gl.
options=odeset('RelTol',10^(-7),'AbsTol',10^(-7));
[t,y]= ode45(Fode,[0,6.65],z0,options);
figure(3);
clf;
plot(y(:,1),y(:,2));
print('-f3','bild3','-dpng','-r100');
%iR Index Loesungsvektor fuer den Abstand Erde kleiner Erdradius
i=1;
while dE(y(i,1),y(i,2))>=RE
        ir=i;
        i=i+1;
end
iR=i
fprintf('Der Flug von der Erde zum Mond und wieder zurueck dauert %6.4f Tage\n',t(iR));
%% Teilaufgabe e
t=t(1:iR,:);
y=y(1:iR,:);
%xKoordinate gegen Zeit
figure(4);
clf;
plot(t,y(:,1));
title('x-Koordinate gegen Zeit');
xlabel('Zeit');
ylabel('x-Koordinate');
print('-f4','bild4','-dpng','-r100');
%Geschwindigkeit in xRichtung gegen Zeit
figure(5);
clf
plot(t,y(:,3));
title('Geschwindigkeit in x-Richtung gegen Zeit');
xlabel('Zeit');
ylabel('Geschwindigkeit in x-Richtung');
print('-f5','bild5','-dpng','-r100');
%Phasendiagramm yKoordinate gegen xKoordinate
figure(6);
clf;
plot(y(:,1),y(:,2));
title('Phasendiagramm');
xlabel('x');
ylabel('y');
print('-f6','bild6','-dpng','-r100');
%Zeitlichen Verlauf des Abstandes der Rakete zur Erde
figure(7);
clf;
plot(t,dE(y(:,1),y(:,2)));
title('Abstand Rakete zur Erde');
xlabel('Zeit');
ylabel('Abstand');
print('-f7','bild7','-dpng','-r100');
%% Teilaufgabe f
%Potential der Rakete als rote Punkte
figure(1); 
hold on;
plot3(y(:,1),y(:,2),pot(y(:,1),y(:,2)),'r.');
hold off;
print('-f1','bild8','-dpng','-r100');
%% Teilaufgabe g
%Minimaler Abstand 
Min=min(dM(y(:,1),y(:,2)))-RM;
minref=Min*3.844*10^8;
fprintf('Minimaler Abstand Rakete Mond %i m\n',minref);
