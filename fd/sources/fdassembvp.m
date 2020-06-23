function [A,b,E,ur] = fdassembvp(grd,opr,bcs,pde)
%% CALL: [A,b,E,ur] = fdassembvp(grd,opr,bcs,pde);
% INPUT:
%    grd ... STRUCT; Gitter-Struktur.
%    opr ... STRUCT; Operatorinformation.
%    bcs ... STRUCT; Randinformationen.
%    pde ... STRUCT; Koeffizientenfunktionen.
% OUTPUT:
%    A ... DOUBLE**; Finite Differenzen Matrix des Problems.
%    b ... DOUBLE*;  Vektor der rechten Seite.
%    E ... DOUBLE**; E ist die Fortsetzung der inneren auf alle Knoten.
%    ur ... DOUBLE*; Vektor mit Dirichlet-Randwerten.
% DESCRIPTION:
% FEMASSEMBVP Assembliert das Randwertproblem.
% [A,B,E,UR] = FEMSOLBVP(GRD,OPR,BCS,PDE) assembliert die reduzierte Matrix A
% und den rechte-Seite-Vektor B des Randwertproblems welches durch die
% Strukturen PDE und BCS durch Anwendung des in OPR definierten
% Finite-Differenzen Verfahrens auf dem Gitter GRD entsteht.
% Um das Problem zu loesen, verwende man etwa V = A\B. V ist dann der reduzierte
% Koeffizientenvektor der Loesung. Mit U = E * V + UR kann man die essentiellen
% Randbedingungen der Loesung zufuegen.
%
% [A,B,E,UR] = FDSOLBVP(GRD,OPR) assembliert das Gleichungssystem des
% Laplace-Randwertproblems mit homogenen Dirichlet Randdaten.
%
% [A,B,E,UR] = FDSOLBVP(GRD,OPR,BCS) assembliert das Gleichungssystem des
% Laplace-Randwertproblems mit Dirichlet Randdaten wie in BCS.

% Version 1.0: Willy Doerfler, KIT, 2013.

% Assemble boundary conditions.
[E,ur] = fdassembcs(grd,opr,bcs);

% Assemble PDE.
[L,f] = fdassempde(grd,opr,pde);

% Reduce the linear system.
A = E'*L*E;
b = E'*(f-L*ur);

return
