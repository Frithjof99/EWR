

clc;
clear;
close all;

loeser = "ode23";

% Definiere Parameter

d = 2.9 * 10 ^ (-2);
a = 2.941 * 10 ^ (-3);
C = 1 - d/(a*3.03);
y0 = 3.03;
t0 = 1960;
t1 = 2200;

yex = @(t) (d / a)./(1 - C * exp(-d*(t - t0)));

verhulst = @(tt,yy) yy * (d - a * yy);

n_cycles = 6;

NN = zeros(n_cycles,1);
err = zeros(n_cycles,1);
ctm = zeros(n_cycles,1);
dt = zeros(n_cycles,1);
figure(1);
clf;
for cycle = 1:n_cycles
   
    opt = odeset('RelTol', 10^(-cycle-1), 'AbsTol', 10^(-cycle-3));
    switch loeser
        case "ode45"
            [t, y] = ode45(verhulst, [t0, t1], y0, opt);
        case "ode23"
            [t, y] = ode23(verhulst, [t0, t1], y0, opt);
    end
    figure(1);
    plot(t,y);
    [a,b] = size(t);
    NN(cycle) = a;
    err(cycle) = max(abs(y-yex(t)));
    
end

figure(2);
clf;
[eoc,cst]=eoctool(NN,err);
