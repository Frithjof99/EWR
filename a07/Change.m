function [dS, dI, dR] = Change(beta, gamma, phi, N_tot, S, I, R)

dS = zeros(27,1);
dR = zeros(27,1);

for i = 1:27
    for j = 1:27
        for k = 1:27
            dS(i) = dS(i) + phi(i,j)*S(i)*(phi(k,j)*I(k)) / N_tot(j);
        end
    end
end
dS = dS * -beta;

for i = 1:27
   dR(i) = gamma * I(i); 
end

dI = - dS - dR;

end