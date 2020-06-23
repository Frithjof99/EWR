% Skript: Run over several number of unknowns
Lierr = zeros(M,1); tm = zeros(M,1);
for jj=1:M
   fprintf('-Run %d of %d\n',jj,M);
   NperDim = N(jj);
   ctm = tic();
   eval(prog);% Run a predefined code
   tm(jj) = max(1e-8,toc(ctm));
   Lierr(jj) = errli;
end
figure(3);
[eoc,cst] = eoctool(N,Lierr);% Experimental order of convergence
aveoc = -mean(eoc(2:end,1)); avcst = mean(cst(2:end,1));
fprintf('\n Geschaetzter Fehler: E = %5.2e * h^{%5.2f}\n',avcst,-aveoc);
figure(4);
[eoc,cst] = eoctool(N,tm);% Ermittlung der Rechenzeit
aveoc = -mean(eoc(2:end,1)); avcst = mean(cst(2:end,1));
fprintf('\n Rechenzeit: T = %5.2e * N^{%5.2f}\n',avcst,aveoc);
title('Growth of cputime');
