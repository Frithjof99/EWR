function [u,fu,ef,out] = newton(f,fp,u0,opts,varargin)
%% CALL: [u,fu,ef,out] = newton(f,fp,u0,opts,varargin);
% INPUT:
%    f,fp ... FUNC; Pointer to f(.,varargin): IR^N->IR^N, resp.
%                   fp(.,varargin): IR^N->IR^{N,N}. By option fpSolve=1, we
%                   use fp(u,fu,varargin)=f'(u)\fu.
%    u0 ... DBLE; Starting value.
%    opts ... STRUCT; Controls
%       opts.Display ... Infolevel: 0=none, 1=initial/end, 2=comments
%                                   3=iteration [DEF 0].
%       opts.MaxIter ... Max. number of iterations [DEF 10].
%       opts.qMax ... Max. number of inner trial steps [DEF 10].
%       opts.StepCtrl ... 0=off, 1=on [DEF 0].
%       opts.TolFun ... f-tolerance [DEF 1.0e-10].
%       opts.TolX ... x-tolerance [DEF 1.0e-8].
%       opts.fpSolve ... Provide mapping f'(u)\fu [Def 0], see above.
%    varargin ... Parameterlist for f,fp, u|->f(u,...).
% OUTPUT:
%    u ... DBLE; Final iterate.
%    fu ... DBLE; f(u,varargin).
%    ef ... INT; Termination code as 'fzero'.
%    out ... STRUCT; Termination structure as for 'fzero'.
% DESCRIPTION:
% Aim: Solve
%    f(u) = 0
% by Newton's method, plain or with Armijo step size control (StepCtrl 0/1).
%%

% Version 2.0. Willy Doerfler, KIT, Mar 29, 2014.

%% Todo:
% Failures, error codes, fp more general

%% Analyse input
switch nargin
case 0% Return default settings
   u.Display = 0;
   u.MaxIter = 10;
   u.qMax = 10;
   u.StepCtrl = 0;
   u.TolFun = 1.0e-12;
   u.TolX = 1.0e-10;
   u.fpSolve = 0;
   return;
case{1,2,3}% Use default settings
   info = 0;
   itmax = 10;
   qmax = 10;
   stepctrl = 0;
   TOLf = 1.0e-12;
   TOLx = 1.0e-10;
   fpSolve = 0;
case {4,5}% User defined settings
   info = opts.Display;
   itmax = opts.MaxIter;
   qmax = opts.qMax;
   stepctrl = opts.StepCtrl;
   TOLf = opts.TolFun;
   TOLx = opts.TolX;
   fpSolve = opts.fpSolve;
end

%% Init
u = u0;
fu = f(u0,varargin{:}); fct = 1;
qq = 0;
if info>=1
   fprintf(' Newton iteration\n');
end

%% Find a zero by (stepsize controlled) Newton's method
for jj=1:itmax
   u_old = u;
   if fpSolve==0% We receive the jacobian
      fpu = fp(u,varargin{:}); fct = fct+1;
      newcorr = fpu\fu;
   else         % We do fp\fu in a subroutine
      newcorr = fp(u,fu,varargin{:}); fct = fct+1;
   end
   if stepctrl==0% No stepsize control
      u = u-newcorr;
      fu = f(u,varargin{:}); fct = fct+1;
   else% Stepsize control (Armijo's rule)
      qq = max(0,qq-1);
      while qq<qmax
         t = 2^(-qq);
         ut = u-t*newcorr;
         fut = f(ut,varargin{:}); fct = fct+1;
         if norm(fut,2)<(1-t/2)*norm(fu,2), break; end
         qq = qq+1;
      end
      u = ut; fu = fut;
   end
   erru = norm(u-u_old,2);
   errf = norm(fu,2);
   if info>=3
      if stepctrl==0
         fprintf('      Nstep %d, errf = %5.2e, erru = %5.2e\n',jj,errf,erru);
      else
         fprintf('      Nstep %2d, errf = %5.2e, erru = %5.2e, qq = %2d\n', ...
                  jj,errf,erru,qq);
      end
   end
   if erru<TOLx || errf<TOLf, break; end
end
if info>=1
   fprintf('   ITstep %d, errf = %7.4e, erru = %7.4e\n',jj,errf,erru);
end

%% Finalize
ef = 1;
out.errx = erru;
out.errf = errf;
out.iterations = jj;
out.funcCount = fct;
if jj>itmax-1
   ef = -7;
   out.message = 'Maximal iteration count reached';
end

%%
end
