%% Federated Streaming Adaptive PCA - Federated Tests
%
% Description: 
%   This code is supplied as additional material alongside our paper:
%   "Federated PCA with Adaptive Rank Estimation"
%   
% This script constructs a virtual federated hierarchy where the PCA 
% object is searched. The result is propageted from bottom to top and the
% final result is stored over at the root node. The result can then be used 
% to globally or selectively enhance modes located at other places by using
% our weighted subspace merging algorithm.
%
% Note:
%  i) Please ensure you have an up-to-date MATLAB version (> 2017a) as 
%     older versions have a problems handling character/string arrays in  
%     certain cases which are extensively used in this script.
% ii) This script uses parallel pool and results may differ depending on
%     CPU & RAM resources available.
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 17/07/2019
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
global allow_print

%% Setup the execution environment

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
% print console flag
allow_print = 0;

% destroy pool at the end of execution
destroy_pool = 0;

% target rank (or seed for adaptive)
r_seed = 10;
% number of features in each vector
feats = 1000; % 1k
% pick the spectrum range to test against
alphas = [0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 2, 3, 4];
%alphas = [0.01, 1, 2];
% the length is the range for testing
alphas_len = size(alphas, 2);
% enable analytical error calculation
no_err_flag = 1;

% Distributed parameters

% desired max tree depths
tree_depths = 1:5;
% the max tree depth
max_tree_depth = tree_depths(end);
% number of feature vectors to process
T_base = (2^max_tree_depth) * 1000; % 2^(max tree depth) * 1k
% t-max
tmax = 4;

% SPCA parameters
fm_rank = r_seed;
fm_blk = 50;
fm_floor_mul = 2;
fm_no_err = no_err_flag;

% array pre-calc
t_total = zeros(tmax, max_tree_depth, alphas_len);
ta_total = zeros(tmax, max_tree_depth, alphas_len);
W_chunk = 0;

% initialise parallel pool
p = gcp;

%% Trial execution

glob_tic = tic; 

% run for the given problem sizes
for t = 1:tmax
  % multiply the problem size by a factor of t
  T = T_base*t;
  fprintf("\n -- Running for T = %d\n", T);
  % tree depth run
  for tree_depth = 1:max_tree_depth
      % set the number chunks (pca workers) for this trial
      chunks = 2^tree_depth;
      % set the chunk size
      chunkSize = T/chunks;
      % clear W_chunk
      clear W_chunk U_chunk G_chunk;
      fprintf("\n -- Running trials for Tree Depth %d", tree_depth);
      for c_alpha = 1:alphas_len
        % run the trials
        p_alpha = alphas(c_alpha);
        fprintf("\n -- Running distributed trial for alpha %d\n",  p_alpha);
        % get the global trial tic
        g_ts = tic;

        % first, create the worker chunks
        c_ts = tic;
        % populate each chunk
        fprintf(['\n ** Generating %d chunks of synthetic data of ', ...
          '%d (feats) x %d (cols)\n'], chunks, feats, chunkSize);
        parfor w = 1:chunks
           W_chunk(w, :, :) = synthetic_data_gen(feats, chunkSize, 1, p_alpha);
        end
        fprintf("\n ** Finished synthetic data generation\n")
        my_toc(c_ts);

        % now for each "global" block compute the chunks and merge
        % them as they become available (in parallel)
        %fprintf(" ** Starting the Federated chunk computation\n");
        svd_ts = tic;
        parfor w = 1:chunks
          [~, ~, U_chunk(w, :, :), G_chunk(w, :, :), ~, ~, ~] = ...
          spca_edge(squeeze(W_chunk(w, :, :)), ...
            fm_rank, fm_blk, ...
            fm_floor_mul, fm_no_err);
        end
        %fprintf(" ** Finished the Federated chunk computation\n");
        pca_tick = my_toc(svd_ts);

        % run for the required levels of merging (which can be done in parallel)
        U_t = U_chunk;
        G_t = G_chunk;
        %fprintf("\n ** Starting Tree Merges\n");
        U_m = NaN(chunks / 2, feats, fm_rank);
        G_m = NaN(chunks / 2, fm_rank, fm_rank);
        for i = 1:tree_depth
          chunk_merges = chunks / (2^i);
          %fprintf("\n  ** STARTED Tree Level %d, requiring %d merges \n", ...
          %  i, chunk_merges);
          level_tic = tic;
          for w = 1:chunk_merges
            %fprintf("\n  -- Worker merging subspaces: %d and %d\n", 2*w-1, 2*w);
            [U_m(w, :, :), G_m(w, :, :)] = ...
              merge_subspaces(squeeze(U_t(2*w-1, :, :)), ...
                              squeeze(G_t(2*w-1, :, :)), ...
                              squeeze(U_t(2*w, :, :)), ...
                              squeeze(G_t(2*w, :, :)));
          end
          %fprintf("\n  ** FINISHED Tree level %d merges \n ", i);
          % swap the variables
          U_t = U_m(1:chunk_merges, :, :);
          G_t = G_m(1:chunk_merges, :, :);
          % the toc level
          merge_tick = my_toc(level_tic);
        end
        %fprintf("\n ** Finished Tree Merges\n");

        % final values for the subspace
        U_final = squeeze(U_t);
        % and the singular values
        G_final = squeeze(G_t);
        % print a nice message
        fprintf("\n -- Finised running distributed trial for alpha %d\n",  ...
          p_alpha);
        % also the time it took to run it
        %t_total(tree_depth, c_alpha) = my_toc(g_ts);
        t_total(t, tree_depth, c_alpha) = (pca_tick + merge_tick);
        ta_total(t, tree_depth, c_alpha) = (pca_tick + merge_tick) / chunks;
      end
      fprintf("\n -- Finished running trials for Tree Depth %d", tree_depth);
  end
  fprintf("\n -- Finished running for T=%d\n", T);
end

%% Print info

% worker scaling
figure;
hold on;
% create the cells required
plegs = cell(1, max_tree_depth);
% make the axes
for i = 1:max_tree_depth
    plegs{i} = sprintf("%d", 2^i);
end
% run all the curves and T legends
wlegs = cell(1, tmax);
for i = 1:tmax
  c_total = squeeze(t_total(i, :, :));
  errorbar(mean(c_total, 2), std(c_total, 0, 2), '-+', 'LineWidth', 2);
  wlegs{i} = sprintf("T=%dK", (T*i) / 1000);
end
% title and stuff
title('Total time execution (workers) for given T');
ylabel('time (s)');
xlabel('node count');
xticks(1:max_tree_depth);
xticklabels(plegs);
legend(wlegs);
% enlarge the fonts
set(gca, 'FontSize', 18);

% amortised worker scaling
figure;
hold on;
% create the cells required
plegs = cell(1, max_tree_depth);
% make the axes
for i = 1:max_tree_depth
    plegs{i} = sprintf("%d", 2^i);
end
% run all the curves and T legends
wlegs = cell(1, tmax);
for i = 1:tmax
  c_total = squeeze(ta_total(i, :, :));
  errorbar(mean(c_total, 2), std(c_total, 0, 2), '-+', 'LineWidth', 2);
  wlegs{i} = sprintf("T=%dK", (T*i) / 1000);
end
% title and stuff
title('Amortised per node execution time (workers) for given T');
ylabel('time (s)');
xlabel('node count');
xticks(1:max_tree_depth);
xticklabels(plegs);
legend(wlegs);
% enlarge the fonts
set(gca, 'FontSize', 18);

%% Clean up

g_time = toc(glob_tic);

% print the total execution time
fprintf("\n ** Total Elapsed time was %d\n", g_time);

% delete the parallel pool
if ~isempty(p) && destroy_pool == 1
  delete(p)
end
