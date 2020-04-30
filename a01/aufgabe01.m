%% Aufgabe 01
% In dieser Aufgabe sollen zwei verschiedene Modelle fuer das
% Verhalten der Weltbevoelkerung verglichen werden.
%
% Autor:   Lena,Hilpp ; Jan Frithjof Fleischhammer
% Version: 30.04.2020

%% Reset Workspace, reset Konsole, alle Fenster schliessen
clear;clc;close all;

%% Initialisierung der benoetigten Fenster
openfigure(2,'init');

%% Aufgabenteil b

% Variablen für exponentielles Wachstum
lambda = 2 * 10^(-2);
t0 = 1960;
c1 = 3.03;

% Variablen für Modell von Verhulst
a = 2.941*10^(-3);
d = 2.9*10^(-2);
c2 = 1 - d/(a*3.03);

% Vektor t
t = 1800:10:2200;

% Ausgabe Schätzwert maximals Anzahl Menschen
fprintf('Schaetzwert Maximum (alt): %5.3e Mrd Menschen\n', d/a);

%% Aufgabenteil c

% Exponentielles Model
y1 = @(t)c1 * exp(lambda*(t - t0));

% Model von Verhulst
y2 = @(t) (d/a)./(1-c2*exp(-d*(t-t0)));

% UN Daten von 1950 bis 2020
tdat=(1950:10:2020);
ydat=[2.54,3.03,3.70,4.46,5.33,6.14,6.96,7.79];

% UN Daten geschaetzt von 2030 bis 2100
test=(2030:10:2100);
yest=[8.55,9.20,9.79,10.15,10.46,10.67,10.81,10.88];

%% Aufgabenteil d

% Graph 1 zeichnen
figure(1);
clf;
hold on;
plot(t, y1(t), 'b:', 'linewidth', 2);
plot(t, y2(t), 'k-', 'linewidth', 2);

plot(tdat, ydat, 'ro');
plot(test, yest, 'r*');

axis([1800 2100 0 12]);

xlabel("Jahr");
ylabel("Mrd. Menschen");

title("Weltbvölkerung 1800 bis 2100");

legend("Exp. Wachstum", "Verhulst", "UN Daten", "UN Schätzung", "location", "SouthEast");

hold off;

print("-f1","bild1","-dpng","-r300");

%% Aufgabenteil e

% Definiere Startwert
y0 = 3.03;

% Definiere parametrisiertes Verhulst-Model
p=[d/a,d];
y3 = @(p,tt) p(1)./(1-(1-p(1)/y0)*exp(-p(2)*(tt-t0)));

% Definiere Least-Square-Fitting Funktion
J = @(p) norm(y3(p,tdat)-ydat,2);

% Startwert für die Problemlösung
p0=[d/a,d];

% Suche Minimum der Least-Square-Fitting Function
pa = fminsearch(J, p0);

% Werte Ausgeben
% J(pa) ist der minimiere, quadratische Abstand zwischen dem
% parametrisierten Verhulst-Model mit berechnetem Parameter und gegebenen
% UN-Werten. Wird klein erwartet
fprintf('J(pa):%6.4f \n',J(pa));
fprintf('Schaetzwert maximale Weltbevoelkerung(neu): %6.4f Mrd. Menschen \n',pa(1));


%% Aufgabenteil f

% Graph 2 zeichnen
figure(2);
clf;
hold on;
plot(t, y2(t), 'k-', 'linewidth', 2);
plot(t, y3(pa,t), 'k:', 'linewidth', 2);

plot(tdat, ydat, 'ro');
plot(test, yest, 'r*');

axis([1800 2100 0 12]);

xlabel("Jahr");
ylabel("Mrd. Menschen");

title("Weltbvölkerung 1800 bis 2100");

legend("Ursprüngliches Verhulst", "Alternatives Verhulst", "UN Daten", "UN Schätzung", "location", "SouthEast");

hold off;

print("-f2","bild2","-dpng","-r300");

%% Aufgabenteil g

% Definiere Funktion, der Nullstelle der Zeitpunkt ist, an dem das Model
% die 8 Mrd Marke überschreitet
y4 = @(t)y3(pa,t)-8;

% Nullstelle Berechnen
ns = fzero(y4, 2010);

% Zeitpunkt in Jahr und Monat unterteilen und ausgeben
year = floor(ns);
month = floor((ns - year) * 12);

fprintf("Es gab etwa ab dem %iten Monat des Jahren %i mehr als 8 Mrd Menschen.", month, year);


