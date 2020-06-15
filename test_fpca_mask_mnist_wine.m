% The script is designed to evaluate the information preservation of F-PCA 
% with and without perturbation masks against the famous MNIST and WINE 
% datasets. Of particular note is the information preservation under 
% perturbation masks as privacy comes at the cost of accuracy. We evaluate 
% across a range of epsilon and delta values which quantatively compare 
% the result against the offline, and F-PCA with no mask applied.
%
% In all instances, for fair comparison, the adaptive rank estimation is 
% disabled for F-PCA.
%
% Based on work of Grammenos et al.: https://arxiv.org/abs/1907.08059
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 02/06/2020
% 
% License: GPLv3
%

%% Initialisation

close all; clear; clc;

% the configuration parameters
rank = 6;
fpca_adaptive = 0;
fpca_private = 1;
fpca_blk_size = 25;
fpca_blk_size_no_mask = 25;

% the dp ranges for epsilon and delta
%epsilons = [0.1, .4];
%deltas = [0.1, 1];

epsilons = [0.1]; % 0.1, 0.6, 1, 2
deltas = [0.1];


% find the length of parameter arrays
epsilons_len = size(epsilons, 2);
deltas_len = size(deltas, 2);

% compute the number of plots
plot_num = epsilons_len * deltas_len;

% dat_desc = "Wine";
dat_desc = "MNIST";

type = sprintf("fpca-mask-%s", dat_desc);

% set the real dataset
params.type = type;
% enable printing
params.pflag = 1;
% use block error
params.use_blk_err = 1;
% print pdfs
params.pdf_print = 1;

% setup the environment
params = setup_vars(params);

if dat_desc == "MNIST"
  dat = csvread('datasets\mnist.csv', 1)';
elseif dat_desc == "Wine"
  dat = csvread('datasets\wine.csv', 1)';
end

% get the colors for the classes
proj_colors = dat(1, :);

% raise if we want to compute the offline PCA
compute_offline = 1;
% raise if we want to compute the MOD-SuLQ
compute_mod_sulq = 1;

%% Offline PCA

% takes a bit of time, skip if we can.
if compute_offline == 1
  % get the left principal subspace
  [Up, ~, ~] = svd(dat);

  % project
  Uproj = Up' * dat;

  % grab the first PC
  first_pc = Uproj(1, :);

  % and the second PC
  second_pc = Uproj(2, :);

  fig = figure;
  % now plot them.
  gscatter(first_pc, second_pc, proj_colors);
  legend('hide');
  
  % titles
  tt = sprintf("Offline PCA %s", dat_desc);
  title(tt);
  ylabel('2nd PC');
  xlabel('1st PC');
  st = sprintf("offline_pca_%s", dat_desc);
  print_fig(fig, st, params);
end

%% F-PCA (no mask)
% set the parameters
fpca_params.adaptive = fpca_adaptive;
fpca_params.blk_size = fpca_blk_size_no_mask;
% run the F-PCA in only adaptive mode (without masks)
[U_fpca, ~, ~] = fpca_edge(dat, rank, fpca_params);

% project
Uproj = U_fpca' * dat;

% grab the first PC
first_pc = Uproj(1, :);

% and the second PC
second_pc = Uproj(2, :);

fig = figure;
% now plot them.
gscatter(first_pc, second_pc, proj_colors);
legend('hide');

% titles
tt = sprintf("F-PCA %s (no mask)", dat_desc);
title(tt);
ylabel('2nd PC');
xlabel('1st PC');

% print the figure
st = sprintf("fpca_no_mask_r_%d_%s", rank, dat_desc);
print_fig(fig, st, params);

%% F-PCA (with perturbation mask)

% handle to the global figure, if needed
fig = figure;
% run the loop to see what we get.
for i = 1:epsilons_len
  for j = 1:deltas_len
    % dial in the current parameters
    fpca_params.private = fpca_private;
    fpca_params.adaptive = fpca_adaptive;
    fpca_params.blk_size = fpca_blk_size;
    fpca_params.holdoff = 5;
    fpca_params.e_p = epsilons(i);
    fpca_params.delta = deltas(j);
    
    % run the fpca in private, adaptive mode
    [U_pfpca, ~] = fpca_edge(dat, rank, fpca_params);

    % project
    Uproj = U_pfpca' * dat;

    % grab the first PC
    first_pc = Uproj(1, :);

    % and the second PC
    second_pc = Uproj(2, :);
    
    % pick the correct subplot in the grid
    subplot(epsilons_len, deltas_len, (i-1) * deltas_len + j);
    % now plot them.
    gscatter(first_pc, second_pc, proj_colors);
    legend('hide');
    
    % titles
    cap = sprintf('F-PCA %s (\\epsilon: %1.2f, \\delta: %1.2f)', ...
      dat_desc, epsilons(i), deltas(j));
    title(cap);
    ylabel('2nd PC');
    xlabel('1st PC');
  end
end

% print the figure
st = sprintf("fpca_with_mask_r_%d_%s", rank, dat_desc);
print_fig(fig, st, params);

%% MOD-SuLQ 
%
% F-PCA should be a bit worse, but not by much.

if compute_mod_sulq == 1
  % handle to the global figure, if needed
  fig = figure;
  % run the loop to see what we get.
  for i = 1:epsilons_len
    for j = 1:deltas_len

      % run mod-sulq
      [U_ms, ~, ~] = svds(mod_sulq(dat, epsilons(i), deltas(j)), rank);

      % project
      Uproj = U_ms' * dat;

      % grab the first PC
      first_pc = Uproj(1, :);

      % and the second PC
      second_pc = Uproj(2, :);

      % pick the correct subplot in the grid
      subplot(epsilons_len, deltas_len, (i-1) * deltas_len + j);
      % now plot them.
      gscatter(first_pc, second_pc, proj_colors);
      legend('hide');
      
      % titles
      cap = sprintf('MOD-SuLQ %s (\\epsilon: %1.2f, \\delta: %1.2f)', ...
        dat_desc, epsilons(i), deltas(j));
      title(cap);
      ylabel('2nd PC');
      xlabel('1st PC');
    end
  end
end

% print the figure
st = sprintf("mod_sulq_r_%d_%s", rank, dat_desc);
print_fig(fig, st, params);