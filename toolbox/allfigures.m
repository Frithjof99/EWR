function of = allfigures(task)
%% Call: of = allfigures(task);
% Clears the set of active figures starting at 1 and that are consecutively
% numbered. It stops when the first non-existing figure is called. There is no
% effect when no figure is open. After activating a figure, the task is
% evaluated. This task may be 'clf'.
% The number 'of' of open figures found is returned.

%% Analyse input
switch nargin
case 0
   task = '';
case 1
otherwise
   error('*** Check input paramter ***');
end
of = 0;

%% Operate on existing figures
for fg=1:10
   try% Existing figures
      get(fg,'number');% Does figure 'fg' exits?
      figure(fg); eval(task); of = of+1;
   catch% Stop if it does not exist
      break;
   end
end
%fprintf(' %d open figures found\n',of);
