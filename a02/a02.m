clear;
clc;
close all;

openfigure(3, 'init');

%% Aufgabe 2

%% Aufgabenteil b

% Definiere Parameter

d = 2.9 * 10 ^ (-2);
a = 2.941 * 10 ^ (-3);
C = 1 - d/(a*3.03);
y0 = 3.03;
t0 = 1960;
t1 = 2200;

yex = @(t) (d / a)./(1 - C * exp(-d*(t - t0)));

verhulst = @(tt,yy) yy * (d - a * yy);

%% Aufgabenteil c

n_cycles = 6;

NN = zeros(n_cycles, 1);
err = zeros(n_cycles, 3);
ctm = zeros(n_cycles, 3);
dt = zeros(n_cycles, 1);

opts = tb_thetaEuler();

opts.theta = 0;

fprintf("Starte Berechnung\n");

figure(1);
clf;
hold on;



hold off;

for cycle = 1:n_cycles
    
    tau = 10 / (2^(cycle-1));
    
    opts.dt = tau;
    figure(1);
    clf;
    hold on;

    t = 1960:10:2200;
    
    plot(t,yex(t),'k--','linewidth',3);
    
    xlabel('Jahre');
    ylabel('Mrd Menschen');
    
    title(strcat('Lösung mit dem theta-Euler-Verfahren (tau=',num2str(tau),')'), 'FontSize', 15);
    
    hold off;
    
    for theta = 0:0.5:1
        
        opts.theta = theta;
        
        tic;
        
        [t, y] = tb_thetaEuler(verhulst, [t0, t1], y0, opts);
        ctm(cycle, 1 + theta * 2) = max(toc, 1e-8);
        figure(1);
        hold on;

        plot(t,y,'linewidth',1);

        hold off;
        drawnow
        
        
        

        dt(cycle) = tau;
        NN(cycle)=ceil((t1-t0)/dt(cycle));
        err(cycle, 1 + theta * 2) = max(abs(y - yex(t)));
    end
    
    legend('exakt', 'theta = 0', 'theta = 0.5', 'theta = 1','location','SouthEast','FontSize',14);
    
    print('-f1',strcat('bild1-',num2str(cycle)),'-dpng','-r100');

    pause(0);
end



%% Teilaufgabe e
%Konvergenzanalyse mittels eoctool
figure(2);
clf;
%eoctool mit TEX Datei
[eoc,cst]=eoctool(NN,err,2,1,{'expl. Euler','Crank-Nicolson','impli. Euler'},'eoc-tabelle.tex');

%gemittelte Werte für eoc(experimentelle Konvergenzordnung) und
%cst(experimentelle Konstante der Fehlerabschaetzung

aveoc = mean(eoc(2:end,1));
avcst = mean(cst(2:end,1));
fprintf('Fehlerschaetzung explizites Eulerverfahren: E = %5.2e * dt^{%5.2f}\n',avcst,aveoc);

aveoc = mean(eoc(2:end,2));
avcst = mean(cst(2:end,2));
fprintf('Fehlerschaetzung Crank-Nicolson-Verfahren: E = %5.2e * dt^{%5.2f}\n',avcst,aveoc);

aveoc = mean(eoc(2:end,3));
avcst = mean(cst(2:end,3));
fprintf('Fehlerschaetzung implizites Eulerverfahren: E = %5.2e * dt^{%5.2f}\n',avcst,aveoc);

print('-f2','bild2','-dpng','-r100');

%% Teilaufgabe f
%statt Fehler Rechenzeiten uebergeben
figure(3);
[eoc,cst] = eoctool(NN,ctm,2,1,{'Rechenzeit expl. Euler', 'Rechenzeit Crank-Nicolson', 'Rechenzeit impli. Euler'},'eoc-Rechenzeiten.tex');
title("Time convergence history");
legend('location','east')
print('-f3','bild3','-dpng','-r100');
