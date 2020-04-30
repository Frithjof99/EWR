function textline(LineNr,LineText)
%% CALL: textline(LineNr,LineText);
% INPUT:
%    LineNr ... INT; Line number for writing.
%    LineText ... STRING; Textstring to write.
% OUTPUT:
%    Writes on the active window.
% DESCRIPTION:
% Simplifies writing into a text window.
%%

% Author: W. Doerfler, KIT, 2015.

LinePos = 1-(LineNr-1)/10;
text(0.0,LinePos,LineText,'FontSize',11);
return
