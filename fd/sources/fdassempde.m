function [L,f] = fdassempde(grd,opr,pde)
%% CALL: [L,f] = fdassempde(grd,opr,pde);
% INPUT:
%    grd ... STRUCT; Gitter-Struktur.
%    opr ... STRUCT; Operatorinformation.
%    pde ... STRUCT; Koeffizientenfunktionen.
% OUTPUT:
%    L ... DOUBLE**; Finite Differenzen Matrix des Problems.
%    f ... DOUBLE*; Vektor der rechten Seite.
% DESCRIPTION:
% FDASSEMPDE Der FD-Stern auf uniformen Finite-Differenzen-Gittern.
% [A,B] = FDASSEMPDE(GRD,OPR,PDE) gibt die Matrix A und den Vektor B der
% rechten Seite des vollen linearen Gleichungssystems zur�ck. Das 
% Gitter ist durch durch die Struktur GRD mit den Komponenten G, X, [Y,]
% NN, NI definiert (siehe FDGRID), die Funktionspointer der partiellen
% Differentialgleichung befinden sich in der Struktur PDE.

% Version 1: Markus Richter, KIT, 2008.
% Version 2: Willy Doerfler, KIT, 2013.

%% Initialisierungen
dim = grd.dim;
pts = opr.stencil;

%% Verteilen
%% DIM=1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dim==1
   if pts==3
      if strcmp(opr.order1,'central')
         [L,f] = fdassempde1d3pt(grd,pde);
      elseif strcmp(opr.order1,'upwind')
         [L,f] = fdassempde1d3pt_upwind(grd,pde);
      end
   end
%% DIM=2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif dim==2
   if pts==5
      if strcmp(opr.order1,'central')
         [L,f] = fdassempde2d5pt(grd,pde);
      elseif strcmp(opr.order1,'upwind')
         [L,f] = fdassempde2d5pt_upwind(grd,pde);
      end
   elseif pts==9
      [L,f] = fdassempde2d9pt(grd,pde);
   else
      error('*** Not implemented ***');
   end
%% DIM=3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif dim==3
   [L,f] = fdassempde3(grd,opr,pde);
else
   error('*** Not implemented ***');
end

end
