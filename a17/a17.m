clc;
close all;
clear;

openfigure(3,'init');
prog = 'statnonlin_problem';% Programmname fuer run_eoc setzen

data = 'dat_bvp2d_nonlin';

lambda = 1;

eval(data);

fddefaults();

lambdas = 0:0.1:6.7;

maximas = zeros(size(lambdas,2),1);



for i = 1:size(lambdas,2)
   lambda = lambdas(i);
   
   eval(data);
   
    statnonlin_problem();
    if mod(lambda, 0.5) == 0 | lambda == 6.7 
    %print('-f2',['bild2-' num2str(lambda) '.png'],'-dpng','-r800');
    end
    maximas(i) = norm(uh,inf);
    
end

print('-f1','bild1','-dpng','-r800');


figure(3);
clf;
plot(lambdas,maximas);


print('-f3','bild3','-dpng','-r800');

