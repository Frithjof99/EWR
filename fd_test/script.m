clear all; close all;
openfigure(4,'init');
prog = 'stationary_problem';% Programmname fuer run_eoc setzen


datafile = 'data';

allfigures('clf');

fddefaults();

eval(datafile);

stationary_problem();

figure(2); title(['Numerische Loesung Bsp(' num2str(1) ')']);

N = [20;40;80;160;320;640];% Knoten pro Dimension
M = length(N);% Anzahl Durchlaeufe
run_eoc;