clear all; close all;clc;
openfigure(4,'init');
prog = 'eigenvalue_problem';% Programmname fuer run_eoc setzen


datafile = 'dat_evp2d_1';

allfigures('clf');

fddefaults();

eval(datafile);

nb = 3;

N = [5;10;20;30;40;50];% Knoten pro Dimension
M = length(N);% Anzahl Durchlaeufe

err = zeros(6,1);

lam_ex = 2*pi^2;

for i = 1:6
   
    eval(datafile);
    NperDim = N(i);
    eigenvalue_problem();
    
    err(i) = min(abs(d-lam_ex));
    
end

figure(4);
[eoc,cst] = eoctool(N,err);% Experimental order of convergence
aveoc = -mean(eoc(2:end,1)); avcst = mean(cst(2:end,1));
fprintf('\n Geschaetzter Fehler: E = %5.2e * h^{%5.2f}\n',avcst,-aveoc);