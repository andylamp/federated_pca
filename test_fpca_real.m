% This script is responsible for evaluating the peformance amongst real
% datasets; four datasets are from the Berkley Motes which have Light,
% Temperature, Volt, and Humidity measurements as well as MNIST and the
% Wine quality datasets. We evaluate the performance of F-PCA (without
% perturbation masks) against the following methods:
%
% - PM: https://arxiv.org/pdf/1307.0032.pdf
% - GROUSE: https://arxiv.org/pdf/1702.01005.pdf
% - SPIRIT: https://dl.acm.org/doi/10.5555/1083592.1083674
% - FD: https://arxiv.org/abs/1501.01711.pdf
%
% Based on work of Grammenos et al.: https://arxiv.org/abs/1907.08059
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 14/06/2020
% 
% License: GPLv3
%

%% Initialisation
clc; clear; close all;

% for reproducibility
rng(300);

% target rank (or seed for adaptive)
r_seed = 10;

% mnist rank
r_mnist = 5;
% wine rank
r_wine = 5;

% err print
err_print = 1;

% set the real dataset
params.type = "real";
% enable printing
params.pflag = 0;
% use block error (faster!)
params.use_blk_err = 1;
% print pdfs
params.pdf_print = 1;

% what to run 
% NOTE: fpca with rank estimation always runs, regardless of options below.
params.fpca_fixed_run = 0;
params.pm_run = 1;
params.sp_run = 1;
params.gr_run = 1;
params.fd_run = 1;

% what to print
params.fro_print = 1;
params.subspace_err_print = 1;
% production print, for shorter titles which are used in the paper
params.prod_print = 1;

% finally, setup the environment
params = setup_vars(params);

% Trial execution

fprintf("\n -- Evaluating F-PCA against real datasets\n\n");

err_comb = [];
leg_comb = {};
idx = 1;

%% MOTES DATASETS

% Light data
light_data = strcat(params.dataset_path, 'q8calibLight.dat');
[lerr, lleg] = eval_fpca_real(light_data, r_seed, 'Light', params);
err_comb = [err_comb; lerr];
leg_comb{idx} = lleg;
idx = idx + 1;

% Temperature data
temp_data = strcat(params.dataset_path, 'q8calibHumTemp.dat');
[terr, tleg] = eval_fpca_real(temp_data, r_seed, 'Temperature', params);
err_comb = [err_comb; terr];
leg_comb{idx} = tleg;
idx = idx + 1;

% Voltage data
volt_data = strcat(params.dataset_path, 'q8calibVolt.dat');
[verr, vleg] = eval_fpca_real(volt_data, r_seed, 'Volt', params);
err_comb = [err_comb; verr];
leg_comb{idx} = vleg;
idx = idx + 1;

% Humidity data
humidData = strcat(params.dataset_path, 'q8calibHumid.dat');
[herr, hleg] = eval_fpca_real(humidData, r_seed, 'Humidity', params);
err_comb = [err_comb; herr];
leg_comb{idx} = hleg;
idx = idx + 1;

%% MNIST

minst_data = strcat(params.dataset_path, 'mnist.csv');
params.is_csv = 1;
[merr, mleg] = eval_fpca_real(minst_data, r_mnist, 'MNIST', params);
err_comb = [err_comb; merr];
leg_comb{idx} = mleg;
idx = idx + 1;

%% WINE

Wine_data = strcat(params.dataset_path, 'wine.csv');
params.is_csv = 1;
[werr, wleg] = eval_fpca_real(Wine_data, r_wine, 'WINE', params);
err_comb = [err_comb; werr];
leg_comb{idx} = wleg;
% if you want to add another dataset below, uncomment this after you add it.
% idx = idx + 1;

%% Print errors (if enabled)

if err_print == 1
  fig = figure;
  % get the necessary details for looping
  xt = 1:size(leg_comb{1}, 2);
  xt_d = 1:size(err_comb(:, 1), 1);
  % plot them in the same graph
  hold on;
  for i = xt
    % else print
    semilogy(xt_d, log(err_comb(:, i)), 'LineWidth', 2);
  end
  hold off;

  % add the legends (should be the same for all in leg_comb)
  legend(leg_comb{1}, 'Location', 'Best');
  
  ylabel('error (log(rmse))');
  xlabel('datasets');
  xticks(xt_d);
  xticklabels({'Light', 'Temp', 'Volt', 'Humid', 'MNIST', 'WINE'});
  title("Subspace Error")
  % increase the font size
  set(gca, 'FontSize', 12);
  st = sprintf("rmse_subspace_err_real_datasets_fd");
  % print the figure
  print_fig(fig, st, params);
end

fprintf("\n -- Finished F-PCA real dataset evaluation\n\n");


