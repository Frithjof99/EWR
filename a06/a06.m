clc;
clear;
close all;

openfigure(3, 'init');

global N;
N = 16;
T = 100;

positions = zeros(3,N);

d = 4*pi / N;

for i = 1:N/2
   
    positions(:,i) = [sin(d*i),cos(d*i),0];
    positions(:,i+N/2) = [sin(d*i),cos(d*i),1];
    
end

figure(1);
plot3(positions(1,:),positions(2,:),positions(3,:),'o');

y0 = reshape(positions,1,[]);

options=odeset('RelTol',10^(-7),'AbsTol',10^(-7));
tic();
[t,y] = ode45(@hydroforce, [0,T],y0, options);
fprintf("Rechenzeit: %f\n", toc());
fprintf("Anzahl Zeitschritte: %i\n", size(t,1));

diffs = diff(t);
fprintf("Größter Zeitschritt: %f\n", max(diffs));
fprintf("Kleinster Zeitschritt: %f\n", min(diffs));

if size(t,1) > 100
    y = interp1(t,y,0:101/100:100);
    t = (0:101/100:100)';
end


figure(2);
clf;
hold on;
plot(y(:,1),y(:,3), 'o', 'Color', [1 0 0]);
plot(y(:,(N*3/2)+1),y(:,(N*3/2)+3), 'o', 'Color', [0 0 0]);
hold off;


c = zeros(16,3);
for i = 1:8
   c(i,:) = [1 0 0];
end
for i = 9:16
   c(i,:) = [0 0 1];
end

for i = 1:size(t,1)
   
    figure(3);
    clf;
    scatter3(y(i,1:3:end),y(i,2:3:end),y(i,3:3:end),16,c,'filled');
    xlim([-3 3]);
    ylim([-3 3]);
    zlim([-60, 5]);
    %pause(0.1);
    
end





