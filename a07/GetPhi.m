function [ Phi ] = GetPhi( alpha )
%GetPhi Berechnet Kopplungsmatrizen
%   Diese Funktion erstellt eine Koppelungsmatrix, in der jeder
% Ortsteil an seine Nachbarn und sich selbst koppelt. Dabei sind die
% Diagonaleinträge $\alpha$ und die Nebendiagonaleintraege $0$ (kein Nachbar) oder
% $\frac{1-\alpha}{N_{Nachbarn}}$. Damit ist die Zeilensumme immer $1$.

Phi = zeros(27,27);
for i=1:27
    Phi(i,i) = alpha;
end
Kopplung{1}  = [27 21];
Kopplung{2}  = [4 5 7 21];
Kopplung{3}  = [19 24];
Kopplung{4}  = [2 6 7 16];
Kopplung{5}  = [2 10 21];
Kopplung{6}  = [4 16];
Kopplung{7}  = [2 4 16 21];
Kopplung{8}  = [23 26 27];
Kopplung{9}  = [11 12 14 16 17 19 22];
Kopplung{10} = [5 15 21 22 24];
Kopplung{11} = [9 12 22 25];
Kopplung{12} = [9 11 14 23 26];
Kopplung{13} = [17 19];
Kopplung{14} = [9 12 17 20 26];
Kopplung{15} = [10 21 22 23 25];
Kopplung{16} = [4 6 7 21];
Kopplung{17} = [9 13 14 19 20];
Kopplung{18} = [8 21 23 27];
Kopplung{19} = [3 9 13 17 24];
Kopplung{20} = [14 17 26];
Kopplung{21} = [1 2 5 7 10 15 16 18 23 27];
Kopplung{22} = [9 10 11 15 24];
Kopplung{23} = [8 12 15 18 21 25 26];
Kopplung{24} = [3 10 19 22];
Kopplung{25} = [11 15 23];
Kopplung{26} = [8 12 14 20 23];
Kopplung{27} = [1 8 18 21];

Ncount = zeros(27,1);
for i=1:27
    cell = Kopplung{i};
    for j=1:size(Kopplung{i},2)
        if cell(j) > 0
            Ncount(i) = Ncount(i) + 1;
        end
    end
end

for i=1:27
    cell = Kopplung{i};
    for j=1:Ncount(i)
        Phi(i, cell(j)) = (1.0-alpha)/Ncount(i);
    end
end

end

