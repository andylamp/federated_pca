%% Test SVD Time Order Independance 
%
% Description: 
%   This code is a toy bench for testing the SVD time order independence 
%   across different permutations of Y for different across different
%   recovery ranks.
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 17/07/2019
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

% destroy pool at the end of execution
destroy_pool = 0;

% target rank (or seed for adaptive)
r = 5;
% number of features in each vector
feats = 100; % 1k
% number of columns
T = 10000;
% pick the spectrum range to test against
alphas = [0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 2, 3, 4];
% lambda for power law
p_lambda = 1;
% the length is the range for testing
alphas_len = size(alphas, 2);
% enable analytical error calculation
no_err_flag = 1;
% trials for hypothesis
max_trials = 20;

% SPCA parameters
fm_rank = r;
fm_blk = 50;
fm_floor_mul = 2;
fm_no_err = no_err_flag;

% PM Parameters
pm_rank = r;
pm_blk = feats;
pm_floor_mul = 2;
pm_no_err = no_err_flag;

% GROUSE Parameters
grouse_rank = r;
gr_no_err = no_err_flag;

% FD Parameters
fd_rank = r;
fd_no_err = no_err_flag;

% SPIRIT Parameters
sp_holdOffTime = 0;
sp_k0 = r;
sp_lambda = .9;
sp_energy = [0.95, 0.98];
sp_silent = 1;
sp_no_err = no_err_flag;

%% Run trials

for a = 1:size(alphas, 2)

  % generate the synthetic data
  Y = synthetic_data_gen(feats, T, p_lambda, alphas(a));

  % run the svds
  [U_n_svds, G_n_svds, ~] = svds(Y, r);

  % run the normal spca
  [~, ~, U_n, G_n, ~, ~, ~] = ...
    spca_edge(Y, ...
      fm_rank, fm_blk, ...
      fm_floor_mul, fm_no_err);

  % run pm (with svds not qr)
  [~, ~, U_n_pm, ~, ~] = ...
    mitliag_pm(Y, pm_rank, pm_blk, pm_floor_mul, pm_no_err);

  % run GROUSE
  [~, ~, U_n_gr, ~, ~] = my_grouse(Y, grouse_rank, gr_no_err);

  % run FD
  [U_n_fd, ~, ~, ~, ~] = fd(Y', fd_rank, fd_no_err);
  U_n_fd = U_n_fd';
  
  % run SPIRIT
  [U_n_sp,  ~, ~, ~, ~, ~, ~, ~] = ...
    SPIRIT(Y', sp_lambda, sp_energy, sp_k0, sp_holdOffTime, sp_silent, ...
      sp_no_err);

  % run for the specified number of trials
  for t = 1:max_trials
    % generate a random permutation
    shuffle_index = randsample(1:T, T);
    % get the matrix for that respective permutation
    Y_p = Y(:, shuffle_index);
    % run spca
    [~, ~, U_p, G_p, ~, ~, ~] = ...
      spca_edge(Y_p, ...
        fm_rank, fm_blk, ...
        fm_floor_mul, fm_no_err);
    % run pm
    [~, ~, U_pm, ~, ~] = ...
      mitliag_pm(Y_p, pm_rank, pm_blk, pm_floor_mul, pm_no_err);  

    % run FD
    [U_fd, ~, ~, ~, ~] = fd(Y', fd_rank, fd_no_err);
    U_fd = U_fd';
  
    % run GROUSE
    [~, ~, U_gr, ~, ~] = my_grouse(Y_p, grouse_rank, gr_no_err);

    % run SPIRIT
    [U_sp,  ~, ~, ~, ~, ~, ~, ~] = ...
      SPIRIT(Y_p', sp_lambda, sp_energy, sp_k0, sp_holdOffTime, sp_silent, ...
        sp_no_err);

    % calculate the mse errors
    u_mse(t, a) = mse(abs(U_p), abs(U_n));
    g_mse(t, a) = mse(abs(G_p), abs(G_n));
    u_s_mse(t, a) = mse(abs(U_p), abs(U_n_svds));
    g_s_mse(t, a) = mse(abs(G_p), abs(G_n_svds));
    pm_u_mse(t, a) = mse(abs(U_pm), abs(U_n_pm));
    pm_s_u_mse(t, a) = mse(abs(U_pm), abs(U_n_svds));
    gr_u_mse(t, a) = mse(abs(U_gr), abs(U_n_gr));
    gr_s_u_mse(t, a) = mse(abs(U_gr), abs(U_n_svds));
    sp_u_mse(t, a) = mse(abs(U_sp), abs(U_n_sp));
    sp_s_u_mse(t, a) = mse(abs(U_sp), abs(U_n_svds));
    fd_u_mse(t, a) = mse(abs(U_fd), abs(U_n_fd));
    fd_s_u_mse(t, a) = mse(abs(U_fd), abs(U_n_svds));
  end

end

% report the errors

% fprintf("\n ** Number Trials %d", t);
% 
% fprintf("\n -- SPCA MSE for U: %d, G: %d", mean(u_mse), mean(g_mse));
% fprintf("\n -- SPCA vs SVDS MSE for U: %d, G: %d\n", ...
%   mean(u_s_mse), mean(g_s_mse));
% 
% fprintf("\n -- PM MSE for U: %d", mean(pm_u_mse));
% fprintf("\n -- PM vs SVDS MSE for U: %d\n", mean(pm_s_u_mse));
% 
% fprintf("\n -- GROUSE MSE for U: %d", mean(gr_u_mse));
% fprintf("\n -- GROUSE vs SVDS for U: %d\n", mean(gr_s_u_mse));
% 
% fprintf("\n -- SPIRIT MSE for U: %d", mean(sp_u_mse));
% fprintf("\n -- SPIRIT vs SVDS for U: %d\n", mean(sp_s_u_mse));

% find the means for the alphas
m_u_mse = mean(u_mse);
%m_g_mse = mean(g_mse);
m_u_s_mse = mean(u_s_mse); 
%m_g_s_mse = mean(g_s_mse);
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



% plot for all the alphas

figure;
hold on;
xt = 1:size(alphas, 2);
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
legend("spca", "spca_s", "pm", "pm_s", "gr", "gr_s", "sp", "sp_s", ...
    "fd", "fd_s");
xticks(xt);
xticklabels(num2cell(alphas));
xlabel("alphas");
ylabel("errors (mse)");
% enlarge the fonts
set(gca, 'FontSize', 12);

figure;
hold on;
xt = 1:size(alphas, 2);
plot(xt, m_u_mse, "g-+", "LineWidth", 2);
%plot(xt, m_g_mse, "LineWidth", 2);
plot(xt, m_u_s_mse, "LineWidth", 2);
%plot(xt, m_g_s_mse, "LineWidth", 2);
plot(xt, m_pm_u_mse, "LineWidth", 2);
plot(xt, m_pm_s_u_mse, "LineWidth", 2);
plot(xt, m_gr_u_mse, "LineWidth", 2);
plot(xt, m_gr_s_u_mse, "LineWidth", 2);
plot(xt, m_sp_u_mse, "LineWidth", 2);
plot(xt, m_sp_s_u_mse, "LineWidth", 2);
plot(xt, m_fd_u_mse, "LineWidth", 2);
plot(xt, m_fd_s_u_mse, "LineWidth", 2);
hold off;
st = sprintf("Errors for subspaces across a (r=%d, feats=%d, T=%dK)", ...
  r, feats, T/1000);
title(st);
legend("spca", "spca_s", "pm", "pm_s", "gr", "gr_s", "sp", "sp_s", ...
    "fd", "fd_s");
xticks(xt);
xticklabels(num2cell(alphas));
xlabel("alphas");
ylabel("errors (mse)");
% enlarge the fonts
set(gca, 'FontSize', 12);

% end error reporting

