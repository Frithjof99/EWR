%% CALL: fddefaults();
% DESCRIPTION:
% FDDEFAULTS besetzt Defaults fuer Koeffizientenfunktionen, Randbedingungen,
% Graphik.

% Version 3.0: Willy Doerfler, KIT, 2017.

%% Datenfunktionen vorbesetzen (siehe Prototyp)
fun_a = 1;% PDE-Koeff. a
fun_b = [];% PDE-Koeff. b (0 dimensionsabhaengig!)
fun_c = 0;% PDE-Koeff. c
fun_f = 0;% PDE-Quellterm f
fun_m = 0;% PDE-Eigenwert-Koeffizient
fun_dir = 0;% Dirichletrandbedingung
fun_xdir = 1;% true, falls Dirichletrandpunkt
fun_xneu = 0;% true, falls Neumannrandpunkt
fun_uex = [];% Default, falls unbekannt
fun_uext = [];% Default, falls unbekannt
% Grafik-Optionen
gopt = fdplot();% Defaultsetzungen fuer fdplot (diskrete Loesung)
gopt.ylbl = 'u_h(x)';% y-Label
gopt2 = fdplot();% Defaultsetzungen fuer fdplot (exakte Loesung)
gopt2.ylbl = 'u(x)';% y-Label
gopt2.ltyp = 'r--';
