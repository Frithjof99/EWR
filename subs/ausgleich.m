function [w,res] = ausgleich(xv,yv,p)
%% CALL: [w,res] = ausgleich(xv,yv,p);
% Adapting data xv|->yv to a polynomial of degree p.
% w contains polynomial coefficients a_0,a_1,... .
% res gives the residual norm(A*xv-yv)
%%

%% Preprocess
% Data is expected in columns. If not, we transpose.
if size(xv,2)~=1, xv = xv'; end
if size(yv,2)~=1, yv = yv'; end

%% Different methods
switch p
case 1
   A = [ones(length(xv),1),xv]; f = yv;
   w = (A'*A)\(A'*f);% Least square solution
otherwise
   error('*** Kein solches p ***');
end

%% Residual (normalised)
res = norm(A*w-f,2)/norm(f,2);
