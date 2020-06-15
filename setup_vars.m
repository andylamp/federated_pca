function [params] = setup_vars(params)
%SETUP_VARS Function responsible for initialising the correct 
% dataset and graph output paths based on OS type while printing
% our current execution configuration.
%
% Based on work of Grammenos et al.: https://arxiv.org/abs/1907.08059
%
% Author Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date 31/05/2020
% 
% License: GPLv3
% 

if ispc
    fprintf("\n !! Detected Windows PC !!\n");
    params.graph_path = ".\graphs\";
    params.dataset_path = ".\datasets\";
else
    fprintf("\n !! Detected Unix-like PC !!\n");
    params.graph_path = './graphs/';
    params.dataset_path = './datasets/';
end

%% Execution parameters

fprintf("\n !! Execution parameters\n");

% check if have a block error enabled
if ~isfield(params, 'use_blk_err')
  params.use_blk_err = 1;
end

% check if we print figures
if ~isfield(params, 'fig_print')
  params.fig_print = 1;
end

% check if we print pdf's
if ~isfield(params, 'pdf_print')
  params.pdf_print = 1;
end

% check if we print png images
if ~isfield(params, 'png_print')
  params.png_print = 1;
end

if params.use_blk_err == 1
  fprintf("\n\t ** Per block error calculation: ENABLED (VERY GOOD!!)");
else
  fprintf("\n\t ** Per block error calculation: DISABLED (snif!)");
end

% if params.run_synthetic == 1
%   fprintf("\n\t ** Running synthetic dataset tests: ENABLED");
% else
%   fprintf("\n\t ** Running synthetic dataset tests: DISABLED");
% end
% 
% if params.run_real == 1
%   fprintf("\n\t ** Running real dataset tests: ENABLED");
% else
%   fprintf("\n\t ** Running real dataset tests: DISABLED");
% end

% check if we have an execution type
if ~isfield(params, 'type')
  fprintf("\n\t ** WARN: No execution type detected, setting to unknown");
  params.type = "unknown";
end

% spacing
fprintf("\n\n");


%% Figure printing configuration

if params.pflag == 1
  fprintf("\n !! Printing functionality is: ENABLED\n");
  % if we do, make sure the directory is created please note that the name
  % adheres to the ISO-8601 timestamp format tagged with the type of run.
  fp = sprintf("%s%s-%s", params.graph_path, ...
    datestr(now, 'yyyymmddTHHMMSS'), params.type);
  fprintf("\n\t** Trying to create save folder in: %s", fp);
  [s, ~, ~] = mkdir(char(fp));
  if s ~= 1
    fprintf("\n\t!! Error, could not create the folder; disabling figure export\n");
    params.pflag = 0;
  else
    fprintf("\n\t** Folder %s created successfully", fp);
    % update the path with the specific runtime folder
    if ispc
        params.graph_path = strcat(fp, "\");
    else
        params.graph_path = strcat(fp, '/');
    end
    fprintf("\n\t** Target graph dir: %s", params.graph_path);
  end
  
  % check if we are printing a .pdf or .png
  if params.pdf_print == 1
      fprintf("\n\t** Printing output is set to: PDF");
  else
      fprintf("\n\t** Printing output is set to: PNG");
  end
  
  % check if we are printing .fig
  if params.fig_print == 1
      fprintf("\n\t** Printing to .fig format is: ENABLED\n");
  else
      fprintf("\n\t** Printing to .fig format is: DISABLED\n");
  end
else
  fprintf("\n !! Printing functionality is: DISABLED\n");
end

end

