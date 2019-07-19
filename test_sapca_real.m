%% Federated Streaming Adaptive PCA - Real Dataset Tests
%
% Description: 
%   This code is supplied as additional material alongside our paper:
%   "Federated PCA with Adaptive Rank Estimation"
%   
% This script is responsible to execute `real_sapca_eval` for all the real
% datasets. Further, it shows some aggregated plots for all datasets 
% processed. The datasets used in our case are the following:
%
%  - Light
%  - Humidity
%  - Temperature
%  - Volt
%
% and were retrived from here: 
%   https://www.cs.cmu.edu/afs/cs/project/spirit-1/www/data/Motes.zip
%
% Notes:
%  i) Please ensure you have an up-to-date MATLAB version (> 2017a) as 
%     older versions have a problems handling character/string arrays in  
%     certain cases which are extensively used in this script.
% ii) This script uses parallel pool and results may differ depending on
%     CPU & RAM resources available.
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 19/07/2019
% 
% License: 
%  code: GPLv3, author: A. Grammenos 
%  paper: A. Grammenos, R. Mendoza-Smith, C. Mascolo, and J. Crowcroft. 
%         Authors retain their respective copyrights. 
%         
%         Pre-print link: https://arxiv.org/abs/1907.08059
%


%% Initialisation
clc; clear; close all;

% scope in globals
global use_blk_err
global pflag
global fig_print
global pdf_print
global datasetPath
global fd_print

% for reproducibility
rng(300);

% use block error
use_blk_err = 0;
% enable printing
pflag = 1;
% use print figs
fig_print = 1;
% print pdfs
pdf_print = 1;
% remove fd for cleaner errors
fd_print = 0;

% target rank (or seed for adaptive)
r_seed = 10;

% err print
err_print = 1;

% setup the environment
setup_vars;

%% Trial execution

fprintf("\n -- Evaluating real datasets\n\n");

% setup path datasets and load them in memory

lightData = strcat(datasetPath, 'q8calibLight.dat');
[lerr, lleg] = real_sapca_eval(lightData, r_seed, 'Light');

tempData = strcat(datasetPath, 'q8calibHumTemp.dat');
[terr, tleg] = real_sapca_eval(tempData, r_seed, 'Temperature');

voltData = strcat(datasetPath, 'q8calibVolt.dat');
[verr, vleg] = real_sapca_eval(voltData, r_seed, 'Volt');

humidData = strcat(datasetPath, 'q8calibHumid.dat');
[herr, hleg] = real_sapca_eval(humidData, r_seed, 'Humidity');

% print errors
if err_print == 1
  figure;
  xt = 1:size(lleg, 2);
  xt_d = 1:4;
  % combine the errors
  err_comb = [lerr; terr; verr; herr];
  fd_idx = find(strcmp(lleg, 'FD'));
  hold on;
  for i = xt
    % check if we print FD
    if fd_print == 1 && ~isempty(fd_idx) && i == fd_idx
      continue;
    end
    % else print
    plot(xt_d, err_comb(:, i), 'LineWidth', 2);
  end
  hold off;

  % add the legends
  pleg = lleg;
  if fd_print == 1 && ~isempty(fd_idx)
   pleg(fd_idx) = [];
  end
  legend(pleg);
  
  ylabel('error (mse)');
  xlabel('datasets');
  xticks(xt_d);
  xticklabels({'Light', 'Temp', 'Volt', 'Humid'});
  title("Subspace Error")
  % increase the font size
  set(gca, 'FontSize', 12);
  
  figure;
  fd_print = 0;
  xt = 1:size(lleg, 2);
  xt_d = 1:4;
  % combine the errors
  err_comb = [lerr; terr; verr; herr];
  fd_idx = find(strcmp(lleg, 'FD'));
  hold on;
  for i = xt
    % check if we print FD
    if fd_print == 1 && ~isempty(fd_idx) && i == fd_idx
      continue;
    end
    % else print
    plot(xt_d, err_comb(:, i), 'LineWidth', 2);
  end
  hold off;

  % add the legends
  pleg = lleg;
  if fd_print == 1 && ~isempty(fd_idx)
   pleg(fd_idx) = [];
  end
  legend(pleg);
  
  ylabel('error (mse)');
  xlabel('datasets');
  xticks(xt_d);
  xticklabels({'Light', 'Temp', 'Volt', 'Humid'});
  title("Subspace Error")
  % increase the font size
  set(gca, 'FontSize', 12);
end

% to end with a nice console offset
fprintf("\n");


