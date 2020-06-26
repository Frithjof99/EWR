function [eoc,cst] = eoctool(NN,EE,gm,eocfac,legstr,outfile)
%% Call 1: [eoc,cst] = eoctool(NN,EE,gm,eocfac,legstr,outfile);
%% Call 2: [eoc,cst] = eoctool(NN,EE,eocopt);
% INPUT:
%    NN ... (N,1); Number of 'unknowns'.
%    EE ... (N,M); Matrix of 'errors' (arbitrary M; <= 4 reasonable).
%    gm ... INT; 'Gliding mean' (fits a subset of 'gm' subsequent
%           datapoints) [DEF 2].
%    eocfac ... INT; To scale the eoc result (see below) [DEF 1].
%    legstr ... CELL; Legend in graphical output [DEF empty].
%    outfile ... INT/CHAR; TEX-table on screen/in this file.
%                Case INT: file id; case CHAR: file name, '': screen
%                [DEF no output].
%    eocopt ... STRUCT; Fieldnames 'gm', 'eocfac', 'legstr', 'outfile',
%               'k0' with meaning and settings as above. Int 'k0' (0...12)
%               adresses a line- and marker style. Get all DEF by calling
%               eocopt = eoctool();
% OUTPUT:
%    eoc, cst ... DOUBLE(N,M); eoc and constant values.
%    Error constants as in (*) below in a formatted table.
% Graphical Output:
%    NN vs EE(:,k) in 'log10-log10' plot
% DESCRIPTION:
% Analysis of (eg) convergence history. NN is assumed to be the number
% characterising the effort (eg, number of unknowns), EE a list of errors
% of different types (eg, EE(:,1) the exact, EE(:,2) the estimated
% error). eoctool now computes the most probable law of the form:
% (*)    EE(:,k) = C_{k,l} NN^{-eoc_{k,l}}
% for some set l. eoc is usually called experimental order of convergence
% (EOC). Often one likes to present dependence on 'h:=NN^{1/spacedim}'. In
% this case use 'eocfac=spacedim'.

% Author: W. Doerfler, Karlsruhe Institute of Technology.
% Last modified: 10.04.2018, 30.04.2020.
%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collection of some further settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intform = '%4';% Integer-Format for NN (omit 'd')
floatform = '%5.2';% Float-Format for EE (omit 'f' or 'e')
str = {'k-o','k--o','k-.o','k:o','b-*','b--*','b-.*','b:*', ...
       'g-x','g--x','g-.x','g:x'};% Linestyle in plot (for <= 12 items)
lw = 2;% LineWidth in plot
ms = 6;% MarkerSize in plot
fsL = 12;% FontSize Label
fsT = 12;% FontSize Title
fsA = 12;% FontSize Axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyse input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch nargin
case 0
   eoc = struct('gm',2,'eocfac',1,'legstr',[],'outfile',[],'k0',0);
   return
case 2
   gm= 2; eocfac= 1; useleg= 0; out= 0; k0= 0;
case 3
   if isstruct(gm)
      eocfac= gm.eocfac;
      legstr= gm.legstr; if isempty(legstr), useleg= 0; end
      outfile= gm.outfile; if isempty(outfile), out= 0; end
      k0 = gm.k0;
      gmn = gm.gm; clear gm; gm = gmn; clear gmn;
   else
      eocfac= 1; useleg= 0; out= 0; k0= 0;
   end
case 4
   useleg= 0; out= 0; k0= 0;
case 5
   useleg= 1; k0= 0;
   if strcmp(char(legstr),''), useleg= 0; end 
   out= 0;
case 6
   useleg= 1; if strcmp(char(legstr),''), useleg= 0; end
   k0= 0;
   if ischar(outfile)
      if strcmp(outfile,''), out= 1; else out= 2; end
   else
      if isnumeric(outfile), out= 3; end
   end
otherwise
   disp('   CALL eoctool(NN,EE,gm,eocfac,legstr,outfile)');
   disp('   At least ''NN,EE'' required');
   return
end
nn = length(NN); mm = size(EE,2);
gmN = min(max(gm,2),nn);
if gmN~=gm, disp('*** Warning in eoctool.m *** gm changed'); end
gm = gmN;
st = nn-gm+1;
if min(NN)<=0, error('*** Error in eoctool.m *** NN not positive'); end
if min(EE)<=0, error('*** Error in eoctool.m *** EE not positive'); end
if mm==1
   fprintf('\n Analyse %d Dataset N vs E with gliding mean %d\n',mm,gm);
else
   fprintf('\n Analyse %d Datasets N vs E with gliding mean %d\n',mm,gm);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute C's and eoc's and produce formatted output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
itf = strcat(intform,'d');
ftf = strcat(floatform,'e'); %ftfE = strcat(floatform,'e');
formstr = [itf ' \t ' ftf ' \t ' ftf ' \t ' ftf '\n'];
eoc = zeros(nn,mm);
cst = zeros(nn,mm);
if st>0
   fprintf(' Error law: e = cst * n^{-eoc}\t (interpolated e)\n\n');
   for k=1:mm
      if useleg==0, fprintf('Data %d\n',k);
      else fprintf('Data %d : %s\n',k,legstr{k}); end
      fprintf('  k \t   e_k \t\t cst_k \t\t eoc_k\n');
      for ll=1:st
         lt = ll+gm-1;
         [est,cst(lt,k),eoc(lt,k)] ...
            = eoctest(NN(ll:ll+gm-1),EE(ll:ll+gm-1,k));
         eoc(lt,k) = eoc(lt,k)*eocfac;
         fprintf(formstr,ll+gm-1,est,cst(lt,k),eoc(lt,k));
      end
   end
else
   error('*** Error in eoctool.m *** Parameter choice does not make sense');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show graphical output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:mm
   plot(log10(NN),log10(EE(:,k)),char(str(k+k0)), ...
        'LineWidth',lw,'MarkerSize',ms);
   hold on;
end
hold off;
title('Convergence history','FontSize',fsT); set(gca,'FontSize',fsA);
xlabel('log_{10}(N)','FontSize',fsL);
ylabel('log_{10}(E)','FontSize',fsL);
if useleg==1, legend(legstr); legend boxoff, end
axis tight; drawnow;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Produce TEX-output
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ftf = strcat(floatform,'f','\\dec{%2d}');
ftfE = strcat(floatform,'f','\\dec{%2d}');
if (out>0)
  sA= '|r|';
  sB= strcat('\Entry',char(65));
  sC= ['$' itf '$'];
  sD= ['$' itf '$'];
  tA= [' & $' ftf '$ '];% The floating strings
  tB= [' & $' ftfE '$'];
  for k=2:2*mm+1
     sA= strcat(sA,'r|');
     sB= strcat(sB,'&','\Entry',char(65+k-1));
     if (mod(k-1,2))
        sC= strcat(sC,tA);
        sD= strcat(sD,tA);
     else
        sC= strcat(sC,' &');
        sD= strcat(sD,tB);
     end
  end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch out
case 1% Show it on the workspace
  fprintf('\nTEX-output from eoctool (edit \\def''s) %s\n\n',date);
  fprintf('\\def%s{N}%% N input\n',strcat('\Entry',char(65)));
  for k=1:mm
     if useleg==0, sS = ''; else sS = legstr{k}; end
     fprintf('\\def%s{%s}%% E%1d input\n',strcat('\Entry',char(65+2*k-1)),sS,k);
     fprintf('\\def%s{eoc}%%\n',strcat('\Entry',char(65+2*k)));
  end
  fprintf('\\def\\dec#1{{}_{#1}}%%\n');
  fprintf('\n\\begin{tabular}{%s}\n\\hline\n',sA);
  fprintf('%s\\\\\n\\hline\n',sB);
  for ll=1:gm-1
     lz = floor(log10(EE(ll,:))); z = EE(ll,:)./10.^lz;
     fprintf(strcat(sC,'\\\\ \n'),NN(ll),reshape([z;lz],2*mm,1));
  end
  for ll=gm:nn
     ee(1:2:2*mm-1)= EE(ll,:); ee(2:2:2*mm)= eoc(ll,:);
     lz= zeros(1,2*mm); lz(1:2:2*mm-1)= floor(log10(ee(1:2:2*mm-1)));
     z = ee./10.^lz;
     fprintf(strcat(sD,'\\\\ \n'),NN(ll),reshape([z;lz],4*mm,1));
  end
  fprintf('\\hline\n\\end{tabular}\n\n');
case {2,3}% Store it in a file
  if out==2% Given filename
     fl= fopen(outfile,'w');
     fprintf('\nStored TEX-output in file "%s"\n\n',outfile);
  elseif out==3% Given file-id
     fl= outfile; 
     fprintf('\nStored TEX-output in file with id %d\n\n',fl);
  end
  fprintf(fl,'TEX-output from eoctool (edit \\def''s) %s\n\n',date);
  fprintf(fl,'\\def%s{N}%%\n',strcat('\Entry',char(65)));
  for k=1:mm
     if useleg==0, sS = ''; else sS = legstr{k}; end
     fprintf(fl,'\\def%s{%s}%%\n',strcat('\Entry',char(65+2*k-1)),sS);
     fprintf(fl,'\\def%s{eoc}%%\n',strcat('\Entry',char(65+2*k)));
  end
  fprintf(fl,'\\def\\dec#1{{}_{#1}}%%\n');
  fprintf(fl,'\n\\begin{tabular}{%s}\n\\hline\n',sA);
  fprintf(fl,'%s\\\\\n\\hline\n',sB);
  for ll=1:gm-1
     lz = floor(log10(EE(ll,:))); z = EE(ll,:)./10.^lz;
     fprintf(fl,strcat(sC,'\\\\ \n'),NN(ll),reshape([z;lz],2*mm,1));
  end
  for ll=gm:nn
     ee(1:2:2*mm-1)= EE(ll,:); ee(2:2:2*mm)= eoc(ll,:);
     lz= zeros(1,2*mm); lz(1:2:2*mm-1)= floor(log10(ee(1:2:2*mm-1)));
     z = ee./10.^lz;
     fprintf(fl,strcat(sD,'\\\\ \n'),NN(ll),reshape([z;lz],4*mm,1));
  end
  fprintf(fl,'\\hline\n\\end{tabular}\n');
  if out==2, fclose(fl); end
otherwise% do nothing
end
% Cut output (avoid zeros for index 1 to gm-1)
eoc = eoc(gm:end,:);
cst = cst(gm:end,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%=======================================================================
%=======================================================================
%=======================================================================

function [errest,cst,eoc] = eoctest(n,err)
% Call:  [errest,cst,eoc] = eoctest(n,err)
% Input:
%    n ... (*,1); Independent variable.
%    err ... (*,1); Error for chosen 'n'.
% Description:
% Analyse data 'n' and 'err' with respect to a assumed law
%    err(i) = cst*n(i)^{-eoc}
% and returns results 'cst' and 'eoc'.
% Method: Take 'log10' and then least square linear fit of data.

% Author: W. Doerfler, Univ. Karlsruhe.
% Last modified:
% 18.03.2002, W. Doerfler. Matlab 6.1.0450.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res = polyfit(log10(n),log10(err),1);
eoc = res(1);
cst = 10^(res(2));
l   = max(size(n));
errest = cst*n(l)^eoc;
eoc    = -eoc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
