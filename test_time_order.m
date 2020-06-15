% The script is designed to evaluate time order independence, in a sense
% how robust the result is against input permutation. The methods evaluated
% against our F-PCA scheme for this metric are the following:
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
% Last touched date: 31/05/2020
% 
% License: GPLv3
%

%% Initialisation
clc; clear; close all;

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
% print mse in console
print_mse_console = 0;

% target rank (or seed for adaptive)
ranks = [5, 20, 40, 60, 80];
% number of features in each vector
feats = 100; % 1k
% number of columns
T = 10000;
% pick the spectrum range to test against
alphas = [0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 2, 3, 4];

% synthetic dataset parameters
synth_params.spectrum_type = "pl";
synth_params.lambda = 1;

% the length is the range for testing
alphas_len = size(alphas, 2);

% trials for hypothesis
max_trials = 20;

% SPIRIT Parameters
sp_lambda = .9;
sp_energy = [0.95, 0.98];

% preallocate error arrays for speed
u_mse = zeros(max_trials, alphas_len);
g_mse = zeros(max_trials, alphas_len);
u_s_mse = zeros(max_trials, alphas_len);
g_s_mse = zeros(max_trials, alphas_len);
pm_u_mse = zeros(max_trials, alphas_len);
pm_s_u_mse = zeros(max_trials, alphas_len);
gr_u_mse = zeros(max_trials, alphas_len);
gr_s_u_mse = zeros(max_trials, alphas_len);
sp_u_mse = zeros(max_trials, alphas_len);
sp_s_u_mse = zeros(max_trials, alphas_len);
fd_u_mse = zeros(max_trials, alphas_len);
fd_s_u_mse = zeros(max_trials, alphas_len);

% set the real dataset
params.type = "time-order";
% enable printing
params.pflag = 1;
% print pdfs
params.pdf_print = 1;

% setup the environment
params = setup_vars(params);



%% Run the trials

for i = 1:size(ranks, 2)
  
  % set the rank
  r = ranks(i);

  fprintf("\n -- Running for rank %d\n", r);
  
  % ---- Run trials -----
  
  for a = 1:alphas_len
    % set the alpha
    synth_params.alpha = alphas(a);
    % generate the synthetic data
    Y = synthetic_data_gen(feats, T, synth_params);

    % run the svds
    [U_n_svds, G_n_svds, ~] = svds(Y, r);

    % run f-pca
    fpca_off_params.blk_size = 50;
    fpca_off_params.adaptive = 0;
    [U_n, G_n, ~] = fpca_edge(Y, r, fpca_off_params);

    % run pm
    [U_n_pm, ~] = mitliag_pm(Y, r);

    % run GROUSE
    [U_n_gr, ~, ~] = my_grouse(Y, r);

    % run FD
    [U_n_fd, ~] = fd(Y', r);
    U_n_fd = U_n_fd';

    % run SPIRIT
    clear sp_params;
    sp_off_params.adaptive = 0;
    sp_off_params.k0 = r;
    sp_off_params.use_qr = 1;
    [U_n_sp,  ~] = SPIRIT(Y', sp_lambda, sp_energy, sp_off_params);

    % run for the specified number of trials
    for t = 1:max_trials
      % generate a random permutation
      shuffle_index = randsample(1:T, T);
      % get the matrix for that respective permutation
      Y_p = Y(:, shuffle_index);

      % run f-pca
      clear fpca_params;
      fpca_params.blk_size = 50;
      fpca_params.adaptive = 0;
      [U_p, G_p, ~] = fpca_edge(Y, r, fpca_params);

      % run pm
      [U_pm, ~] = mitliag_pm(Y, r);

      % run GROUSE
      [U_gr, ~, ~] = my_grouse(Y, r);

      % run SPIRIT
      clear sp_params;
      sp_params.adaptive = 0;
      sp_params.k0 = r;
      sp_params.use_qr = 1;
      [U_sp,  ~] = SPIRIT(Y', sp_lambda, sp_energy, sp_params);

      % run FD
      [U_fd, ~] = fd(Y', r);
      U_fd = U_fd';

      % calculate the mse errors
      u_mse(t, a) = log(immse(abs(U_p), abs(U_n)));
      g_mse(t, a) = log(immse(abs(G_p), abs(G_n)));
      u_s_mse(t, a) = log(immse(abs(U_p), abs(U_n_svds)));
      g_s_mse(t, a) = log(immse(abs(G_p), abs(G_n_svds)));
      pm_u_mse(t, a) = log(immse(abs(U_pm), abs(U_n_pm)));
      pm_s_u_mse(t, a) = log(immse(abs(U_pm), abs(U_n_svds)));
      gr_u_mse(t, a) = log(immse(abs(U_gr), abs(U_n_gr)));
      gr_s_u_mse(t, a) = log(immse(abs(U_gr), abs(U_n_svds)));
      sp_u_mse(t, a) = log(immse(abs(U_sp), abs(U_n_sp)));
      sp_s_u_mse(t, a) = log(immse(abs(U_sp), abs(U_n_svds)));
      fd_u_mse(t, a) = log(immse(abs(U_fd), abs(U_n_fd)));
      fd_s_u_mse(t, a) = log(immse(abs(U_fd), abs(U_n_svds)));
    end

  end

  % find the means for the alphas
  m_u_mse = mean(u_mse);
  m_g_mse = mean(g_mse);
  m_u_s_mse = mean(u_s_mse); 
  m_g_s_mse = mean(g_s_mse);
  m_pm_u_mse = mean(pm_u_mse);
  m_pm_s_u_mse = mean(pm_s_u_mse);
  m_gr_u_mse = mean(gr_u_mse);
  m_gr_s_u_mse = mean(gr_s_u_mse);
  m_sp_u_mse = mean(sp_u_mse);
  m_sp_s_u_mse = mean(sp_s_u_mse);
  m_fd_u_mse = mean(fd_u_mse);
  m_fd_s_u_mse = mean(fd_s_u_mse);

  % and the std's
  std_u_mse = std(u_mse);
  %std_g_mse = std(g_mse);
  std_u_s_mse = std(u_s_mse); 
  %std_g_s_mse = std(g_s_mse);
  std_pm_u_mse = std(pm_u_mse);
  std_pm_s_u_mse = std(pm_s_u_mse);
  std_gr_u_mse = std(gr_u_mse);
  std_gr_s_u_mse = std(gr_s_u_mse);
  std_sp_u_mse = std(sp_u_mse);
  std_sp_s_u_mse = std(sp_s_u_mse);
  std_fd_u_mse = std(fd_u_mse);
  std_fd_s_u_mse = std(fd_s_u_mse);

  % report the errors

  % check if we print info in console
  if print_mse_console == 1
    fprintf("\n ** Number Trials %d", t);

    fprintf("\n -- F-PCA MSE for U: %d, G: %d", m_u_mse, m_g_mse);
    fprintf("\n -- F-PCA vs SVDS MSE for U: %d, G: %d\n", m_u_s_mse, m_g_s_mse);

    fprintf("\n -- PM MSE for U: %d", m_pm_u_mse);
    fprintf("\n -- PM vs SVDS MSE for U: %d\n", m_pm_s_u_mse);

    fprintf("\n -- GROUSE MSE for U: %d", m_gr_u_mse);
    fprintf("\n -- GROUSE vs SVDS for U: %d\n", m_gr_s_u_mse);

    fprintf("\n -- SPIRIT MSE for U: %d", m_sp_u_mse);
    fprintf("\n -- SPIRIT vs SVDS for U: %d\n", m_sp_s_u_mse);

    fprintf("\n -- FD MSE for U: %d", m_fd_u_mse);
    fprintf("\n -- FD vs SVDS for U: %d\n", m_fd_s_u_mse);
  end


  % ---- Plot the results -----

  fig = figure;
  xt = 1:alphas_len;
  hold on;
    errorbar(xt, m_u_mse, std_u_mse, "g-+", "LineWidth", 2);
    %errorbar(xt, m_g_mse, std_g_mse, "LineWidth", 2);
    errorbar(xt, m_u_s_mse, std_u_s_mse, "LineWidth", 2);
    %errorbar(xt, m_g_s_mse, std_g_s_mse, "LineWidth", 2);
    errorbar(xt, m_pm_u_mse, std_pm_u_mse, "LineWidth", 2);
    errorbar(xt, m_pm_s_u_mse, std_pm_s_u_mse, "LineWidth", 2);
    errorbar(xt, m_gr_u_mse, std_gr_u_mse, "LineWidth", 2);
    errorbar(xt, m_gr_s_u_mse, std_gr_s_u_mse, "LineWidth", 2);
    errorbar(xt, m_sp_u_mse, std_sp_u_mse, "LineWidth", 2);
    errorbar(xt, m_sp_s_u_mse, std_sp_s_u_mse, "LineWidth", 2);
    errorbar(xt, m_fd_u_mse, std_fd_u_mse, "LineWidth", 2);
    errorbar(xt, m_fd_s_u_mse, std_fd_s_u_mse, "LineWidth", 2);
  hold off;
  st = sprintf(['Errors (with bars) for subspaces across a ', ... 
    '(r=%d, feats=%d, T=%dK)'], r, feats, T/1000);
  title(st);
  legend("fpca", "fpca_s", "pm", "pm_s", "gr", "gr_s", "sp", "sp_s", ...
      "fd", "fd_s");
  xticks(xt);
  xticklabels(num2cell(alphas));
  xlabel("\alpha");
  ylabel("errors (log(rmse))");
  
  % generate the specific filename
  st = sprintf("time_order_ind_r_%d_d_%d_T_%d_log_rmse", r, feats, T);
  % print the figure
  print_fig(fig, st, params);
  
  % enlarge the fonts
  set(gca, 'FontSize', 12);
  
  fprintf(" -- Finished running for rank %d", r);

  % end error reporting, move on to the next r
end

