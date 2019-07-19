function [] = print_fig(fig, fname, ptype)
%PRINT_FIG This function is responsible for printing the figure to 
% `fname` using `ptype` (default is `-dpng`)
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 30/12/2018
% 
% License: GPLv3
%

% let the function know of these global vars
global pflag
global pdf_print
global fig_print
global graphPath

% check if we have a ptype argument
if nargin < 3
  ptype = '-dpng'; % change to "-dpdf" for pdf
end

% check if we are allowed to print
if pflag == 1
  % generate the full path by concat the graph path + fname
  p = strcat(graphPath, fname);
  % check if we are saving in a .fig format
  if fig_print == 1
    savefig(fig, char(p));
  end
  % print normally (as an image)
  print(fig, char(p), ptype);
  
  % check if we have a global flag that we print to pdf as well
  if pdf_print == 1
    ptype = '-dpdf'; % change to '-dpdf' due to global flag
    print(fig, char(p), ptype);
  end
end

end

