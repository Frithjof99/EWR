domain = 'I';
geom = [1/2 1];
NperDim = 65;

star = 3;
order1 = 'central';

% Datenfunktionen
   fun_f = 2;
   fun_dir = 0;
   fun_uex = @(x) x.*(1-x);
   % Grafik-Optionen
   umin = 0; umax = 0.3;
   
   gopt.fixaxis = [geom(1)-geom(2)/2 geom(1)+geom(2)/2 umin umax];