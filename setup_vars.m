function [] = setup_vars()
%SETUP_VARS If a function responsible for initialising the correct 
% dataset and graph output paths based on OS type while printing
% our current execution configuration.
%
% Author Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date 06/06/2018
% 
% License: GPLv3
% 

% invoke the global variables
global graphPath
global pflag
global pdf_print
global fig_print
global datasetPath

global run_synthetic
global run_real
global use_blk_err

if ispc
    fprintf("\n !! Detected Windows PC !!\n");
    graphPath = ".\graphs\";
    datasetPath = ".\datasets\";
else
    fprintf("\n !! Detected Unix-like PC !!\n");
    graphPath = './graphs/';
    datasetPath = './datasets/';
end

%% Execution parameters

fprintf("\n !! Execution parameters\n");

if use_blk_err == 1
  fprintf("\n\t ** Per block error calculation: ENABLED (VERY GOOD!!)");
else
  fprintf("\n\t ** Per block error calculation: DISABLED (snif!)");
end

if run_synthetic == 1
  fprintf("\n\t ** Running synthetic dataset tests: ENABLED");
else
  fprintf("\n\t ** Running synthetic dataset tests: DISABLED");
end

if run_real == 1
  fprintf("\n\t ** Running real dataset tests: ENABLED");
else
  fprintf("\n\t ** Running real dataset tests: DISABLED");
end

% spacing
fprintf("\n\n");


%% Figure printing configuration

if pflag == 1
  fprintf("\n !! Printing functionality is: ENABLED\n");
  % if we do, make sure the directory is created please note that the name
  % adheres to the ISO-8601 timestamp format
  fp = strcat(graphPath, datestr(now, 'yyyymmddTHHMMSS'));
  fprintf("\n\t** Trying to create save folder in: %s", fp);
  [s, ~, ~] = mkdir(char(fp));
  if s ~= 1
    fprintf("\n\t!! Error, could not create the folder; disabling figure export\n");
    pflag = 0;
  else
    fprintf("\n\t** Folder %s created successfully", fp);
    % update the path with the specific runtime folder
    if ispc
        graphPath = strcat(fp, "\");
    else
        graphPath = strcat(fp, '/');
    end
    fprintf("\n\t** Target graph dir: %s", graphPath);
  end
  
  % check if we are printing a .pdf or .png
  if pdf_print == 1
      fprintf("\n\t** Printing output is set to: PDF");
  else
      fprintf("\n\t** Printing output is set to: PNG");
  end
  
  % check if we are printing .fig
  if fig_print == 1
      fprintf("\n\t** Printing to .fig format is: ENABLED\n");
  else
      fprintf("\n\t** Printing to .fig format is: DISABLED\n");
  end
else
  fprintf("\n !! Printing functionality is: DISABLED\n");
end

end

