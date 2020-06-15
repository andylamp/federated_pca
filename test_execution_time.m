% The script is designed to evaluate the speed of F-PCA against other
% competing methods - specifically we compare against:
%
% - PM: https://arxiv.org/pdf/1307.0032.pdf
% - GROUSE: https://arxiv.org/pdf/1702.01005.pdf
% - SPIRIT: https://dl.acm.org/doi/10.5555/1083592.1083674
% - FD: https://arxiv.org/abs/1501.01711.pdf
%
% The following parameters can be configured:
%
% - n: the ambient dimension
% - r: r-truncation target
% - alpha: parameter for the power law distribution
% - trials: the number of performance runs for each parameter tuple.
%
% We define the above as a parameter tuple and the defaults are the
% following: (200:200:1000, 10, 1, 5)
%
% Based on work of Grammenos et al.: https://arxiv.org/abs/1907.08059
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 08/06/2020
% 
% License: GPLv3
%

%% Initialisation

% fist clear everything
clear; close all; clc;

% trial parameters
d = 200:200:1000;
r = 10;
trials = 3;
T = 0;
use_log_plot = 0;

% dataset parameters
synth_params.spectrum_type = "pl";
synth_params.lambda = 1;
synth_params.alpha = 1;

% general parameters
params.pflag = 1;
params.type = "speed-test";

% setup the parameters
params = setup_vars(params);

% raise the no error flag for speed tests
no_err = 1;
% number of parameter tuples
param_num = size(d, 2);

% preallocate time arrays
fpca_ta = zeros(param_num, trials);
fpca_mask_ta = zeros(param_num, trials);
pm_ta = zeros(param_num, trials);
gr_ta = zeros(param_num, trials);
fd_ta = zeros(param_num, trials);
sp_ta = zeros(param_num, trials);

if T < 2000
  fprintf("\n ** WARN: T must be at least 2k, reverting to default 10k\n");
  T = 10000;
end

%% Run the trials

% now run the speed test
for i = 1:param_num
  fprintf("\n ** Starting speed %d trials for: d=%d, r=%d\n", trials, d(i), r);
  % run the number of specified trials
  for j = 1:trials
    fprintf("\n -- Trial no: %d...\n", j);
    % first, generate the synthetic data for that trial
    Y = synthetic_data_gen(d(i), T, synth_params);
    
    % secondly, execute the timed runs
    clear params_fpca
    params_fpca.verbose = 0;
    [~, ~, fpca_opts] = fpca_edge(Y, r, params_fpca);
    clear params_fpca_mask
    params_fpca_mask.verbose = 0;
    params_fpca_mask.private = 1;
    [~, ~, fpca_mask_opts] = fpca_edge(Y, r, params_fpca_mask);
    [~, pm_opts] = mitliag_pm(Y, r);
    [~, ~, gr_opts] = my_grouse(Y, r);
    [~, fd_opts] = fd(Y', r);
    clear params_sp
    params_sp.verbose = 0;
    [~, sp_opts] = SPIRIT(Y', 0.95, [0.95,0.98], params_sp);
    
    % then, set the execution time for this trial
    fpca_ta(i, j) = fpca_opts.t;
    fpca_mask_ta(i, j) = fpca_mask_opts.t;
    pm_ta(i, j) = pm_opts.t;
    gr_ta(i, j) = gr_opts.t;
    fd_ta(i, j) = fd_opts.t;
    sp_ta(i, j) = sp_opts.t;
    
    % clear the structures
    clear params_fpca;
    clear params_fpca_mask;
    clear params_sp;
  end
  fprintf("\nFinished running %d trials...\n", j);
end

%% Plot the results

% find the means for the speed runs

if use_log_plot == 1
  mtr_fpca = log(mean(fpca_ta, 2));
  mtr_mask_fpca = log(mean(fpca_mask_ta, 2));
  mtr_pm = log(mean(pm_ta, 2));
  mtr_gr = log(mean(gr_ta, 2));
  mtr_fd = log(mean(fd_ta, 2));
  mtr_sp = log(mean(sp_ta, 2));

  % find the std for the speed runs
  str_fpca = log(std(fpca_ta, 0, 2));
  str_mask_fpca = log(std(fpca_mask_ta, 0, 2));
  str_pm = log(std(pm_ta, 0, 2));
  str_gr = log(std(gr_ta, 0, 2));
  str_fd = log(std(fd_ta, 0, 2));
  str_sp = log(std(sp_ta, 0, 2));
else
  % find the means for the speed runs
  mtr_fpca = mean(fpca_ta, 2);
  mtr_mask_fpca = mean(fpca_mask_ta, 2);
  mtr_pm = mean(pm_ta, 2);
  mtr_gr = mean(gr_ta, 2);
  mtr_fd = mean(fd_ta, 2);
  mtr_sp = mean(sp_ta, 2);

  % find the std for the speed runs
  str_fpca = std(fpca_ta, 0, 2);
  str_mask_fpca = std(fpca_mask_ta, 0, 2);
  str_pm = std(pm_ta, 0, 2);
  str_gr = std(gr_ta, 0, 2);
  str_fd = std(fd_ta, 0, 2);
  str_sp = std(sp_ta, 0, 2);
end

%% F-PCA (without mask)

% plot the figure
fig = figure;
hold on
  errorbar(mtr_sp, str_fd, '-^', 'LineWidth', 2);
  errorbar(mtr_pm, str_pm, '-*', 'LineWidth', 2);
  errorbar(mtr_fd, str_fd, '-x', 'LineWidth', 2);
  errorbar(mtr_gr, str_gr, '-+', 'LineWidth', 2);
  errorbar(mtr_fpca, str_fpca, '-o', 'LineWidth', 2);
hold off

% full legend cells
legendCells = {'SP', 'PM', 'FD', 'GROUSE', 'F-PCA'}; 

% assign labels
legend(legendCells, 'location', 'best');
xticks(1:1:size(d, 2));
xticklabels(num2cell(d));
xlabel("ambient dimension (d)"); 
if use_log_plot == 1
  ylabel("average per trial time (log(sec))");
else
  ylabel("average per trial time (sec)");
end
cap = sprintf("Time to compute U for r=%d", r);
title(cap);

% print figure, if needed
t = sprintf("speedtest_T_%sk_kr_%d_alpha_%d_trials_%d", ...
  strrep(num2str(T/1000), ".", "_"), r, synth_params.alpha, trials);
print_fig(fig, t, params);

%% Differentially Private F-PCA (with mask) included

% plot the figure
fig = figure;
hold on
  errorbar(mtr_sp, str_sp, '-^', 'LineWidth', 2);
  errorbar(mtr_pm, str_pm, '-*', 'LineWidth', 2);
  errorbar(mtr_fd, str_fd, '-x', 'LineWidth', 2);
  errorbar(mtr_gr, str_gr, '-+', 'LineWidth', 2);
  errorbar(mtr_fpca, str_fpca, '-o', 'LineWidth', 2);
  errorbar(mtr_mask_fpca, str_mask_fpca, '--', 'LineWidth', 2);
hold off

% full legend cells
legendCells = {'SP', 'PM', 'FD', 'GROUSE', 'F-PCA', 'F-PCA_{mask}'}; 

% assign labels
legend(legendCells, 'location', 'best');
xticks(1:1:size(d, 2));
xticklabels(num2cell(d));
xlabel("ambient dimension (d)");  
if use_log_plot == 1
  ylabel("average per trial time (log(sec))");
else
  ylabel("average per trial time (sec)");
end
cap = sprintf("Time to compute U for r=%d", r);
title(cap);

% print figure, if needed
t = sprintf("speedtest_mask_T_%sk_kr_%d_alpha_%d_trials_%d", ...
  strrep(num2str(T/1000), ".", "_"), r, synth_params.alpha, trials);
print_fig(fig, t, params);