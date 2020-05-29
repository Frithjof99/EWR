function F = hydroforce(t,y)

global N;

K = @(z) (eye(3)+z*z'./(norm(z)^2))./norm(z);

e = [0,0,1]';

Y = reshape(y,3,[]);

F = zeros(3*N,1);

for alpha = 1:N
   r = 0;
   for beta = 1:N
      if beta == alpha
          continue;
      end
      
      x_alpha = Y(:,alpha);
      x_beta = Y(:,beta);
                  
      r = r + K(x_alpha-x_beta)*e;
      
   end
   F((alpha-1)*3+1:(alpha-1)*3+3) = r * (-5) ./ (8*N);
end