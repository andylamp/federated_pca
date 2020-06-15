function [] = print_fig(fig, fname, params)
%PRINT_FIG This function is responsible for printing the figure to 
% `fname` using `ptype` (default is `-dpng`)
%
% Based on work of Grammenos et al.: https://arxiv.org/abs/1907.08059
%
% Author Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date 31/05/2020
% 
% License: GPLv3
% 

% check if we have a print flag (default is disabled)
if ~isfield(params, 'pflag')
  pflag = 0;
else
  pflag = params.pflag;
end

% check if we have a pdf print flag (default is enabled)
if ~isfield(params, 'pdf_print')
  pdf_print = 1;
else
  pdf_print = params.pdf_print;
end

% check if we have a fig print
if ~isfield(params, 'fig_print')
  fig_print = 1;
else
  fig_print = params.fig_print;
end

% check if we have a graph path set
if ~isfield(params, 'graph_path')
  error('** ERR: Graph path must be present')
else
  graph_path = params.graph_path;
end

% check if we are allowed to print
if pflag == 1
  % generate the full path by concat the graph path + fname
  p = strcat(graph_path, fname);
  % check if we are saving in a .fig format
  if fig_print == 1
    savefig(fig, char(p));
  end
  % print normally (as an image)
  print(fig, char(p), '-dpng');
  
  % check if we have a global flag that we print to pdf as well
  if pdf_print == 1
    ptype = '-dpdf'; % change to '-dpdf' for pdf prints
    print(fig, char(p), ptype);
  end
end

end

