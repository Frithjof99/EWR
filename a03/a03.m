%% Reset Workspace, Ausgabe organisieren

clear;
clc;
close all;

openfigure(4, 'init');


%% Raeuber-Beute-Modell (Lotka-Volterra)
ar = 2;
ab = 1;
gr = 1;
gb = 1;
r0 = 2;
b0 = 2.25;

ym = @(tt, yy)[-ar*yy(1)+gr*yy(2)*yy(1),ab*yy(2)-gb*yy(1)*yy(2)]';
y0 = [r0, b0];

%% ODE-Loeser
figure(1);
clf;

optEuler = tb_thetaEuler;
optEuler.theta = 1;
optEuler.dt = 0.025;
optEuler.ssolver = 'fminsearch';

for i = 1:2
   
    switch i
        case 1
            fprintf("ode45\n");
            [t, y] = ode45(ym,[0,10], y0);
        case 2
            fprintf("euler\n");
            [t, y] = tb_thetaEuler(ym, [0,10],y0, optEuler);
            y = y';
    end
    
    figure(1);
    hold on;
    
    
    
    plot(t,y(:,1),'r','linewidth', 2);
    plot(t,y(:,2),'b','linewidth', 2);
    
    
    hold off;
    
    fprintf("Minimale Räuber: %f\n", min(y(:,1)));
    fprintf("Minimale Beute: %f\n", min(y(:,2)));
    fprintf("Maximale Räuber: %f\n", max(y(:,1)));
    fprintf("Maximale Beute: %f\n", max(y(:,2)));
    
    figure(2);
    hold on;
    
    plot(y(:,1), y(:,2));
    
    hold off;
    
    %% Phasendiagramm mit Potential

    f = @(rr,bb)exp(ab*log(rr)-gb*rr+ar*log(bb)-gr*bb);
    
    x = 0:0.1:3;
    z = 0:0.1:4;
    
    [X, Y] = meshgrid(x,z);
        
    figure(3);
        
    hold on;
    contour(x,z,f(X,Y));
    plot(y(:,1), y(:,2));
    hold off;
end













