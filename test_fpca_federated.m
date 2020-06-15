% This script is responsible to perform the federated tests; it creates a
% large enough datasets and arranges the nodes based on the selected tree
% depth. Due to matlabs multiprocessing limitations we can only use a
% binary tree structure but this is not required but rather selected for
% convenience.
%
% The end result shows how much processing time it takes to perform the
% actual PCA as well as the respective merges as the answer is propagated
% upwards. Both amortised and actual times are reported for total, pca
% computation, and merging times.
%
%
% Based on work of Grammenos et al.: https://arxiv.org/abs/1907.08059
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 03/06/2020
% 
% License: GPLv3
%

%% Initialisation
clc; clear; close all;

% for reproducibility
rng(300);

% destroy pool at the end of execution
destroy_pool = 0;

% target rank (or seed for adaptive)
rank = 10;
% number of features in each vector
feats = 1000; % 1k

% synthetic dataset parameters
synth_params.spectrum_type = "pl";
synth_params.lambda = 1;

% pick the spectrum range to test against
alphas = [0.01, 0.1, .5, 1];
% the length is the range for testing
alphas_len = size(alphas, 2);

% Federation parameters

% desired max tree depths
tree_depths = 1:6; % 6
% the max tree depth
max_tree_depth = tree_depths(end);
% number of feature vectors to process
T_base = (2^max_tree_depth) * 2000; % 2^(max tree depth) * 1k
% t-max
tmax = 5; % 5

% chunk number
W_chunk = 0;


% F-PCA parameters
fpca_opts.adaptive = 0;

% array pre-calc
t_total = zeros(tmax, max_tree_depth, alphas_len);
ta_total = zeros(tmax, max_tree_depth, alphas_len);

merge_time = zeros(tmax, max_tree_depth, alphas_len);
pca_time = zeros(tmax, max_tree_depth, alphas_len);

merge_time_ta = zeros(tmax, max_tree_depth, alphas_len);
pca_time_ta = zeros(tmax, max_tree_depth, alphas_len);

% setup the variables
params.type = "federated-test";
% enable printing
params.pflag = 1;
% print pdfs
params.pdf_print = 1;

% create the environment parameters
params = setup_vars(params);

% initialise parallel pool
p = gcp;

%% Federated Trial execution

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
        % assign the alpha to the parameters
        synth_params.alpha = p_alpha;
        fprintf("\n -- Running federated trial for alpha %d\n",  p_alpha);
        % get the global trial tic
        g_ts = tic;

        % first, create the worker chunks
        c_ts = tic;
        % populate each chunk
        fprintf(['\n ** Generating %d chunks of synthetic data of ', ...
          '%d (feats) x %d (cols)\n'], chunks, feats, chunkSize);
        parfor w = 1:chunks
           W_chunk(w, :, :) = synthetic_data_gen(feats, chunkSize, synth_params);
        end
        fprintf("\n ** Finished synthetic data generation\n")
        my_toc(c_ts);

        % now for each "global" block compute the chunks and merge
        % them as they become available (in parallel)
        
        % Chunk computation
        fprintf(" ** Starting the Federated chunk computation\n");
        
        svd_ts = tic;
        parfor w = 1:chunks
          [U_chunk(w, :, :), G_chunk(w, :, :), ~] = ...
            fpca_edge(squeeze(W_chunk(w, :, :)), rank, fpca_opts);
        end
        pca_tick = my_toc(svd_ts);
        
        % Finished chunk computation
        fprintf(" ** Finished the Federated chunk computation\n");
        

        % run for the required levels of merging (which can be done in parallel)
        U_t = U_chunk;
        G_t = G_chunk;
        
        % Merge computation
        fprintf("\n ** Starting: Merges for tree level %d\n", i);
        
        U_m = NaN(chunks / 2, feats, rank);
        G_m = NaN(chunks / 2, rank, rank);
        merge_tick = 0;
        for i = 1:tree_depth
          chunk_merges = chunks / (2^i);
          % fprintf("\n  ** STARTED Tree Level %d, requiring %d merges \n", ...
          %  i, chunk_merges);
          level_tic = tic;
          for w = 1:chunk_merges
            % fprintf("\n  -- Worker merging subspaces: %d and %d\n", 2*w-1, 2*w);
            [U_m(w, :, :), G_m(w, :, :)] = ...
              fpca_subspace_merge(squeeze(U_t(2*w-1, :, :)), ...
                                  squeeze(G_t(2*w-1, :, :)), ...
                                  squeeze(U_t(2*w, :, :)), ...
                                  squeeze(G_t(2*w, :, :)));
          end
          % add the time for the current level
          m_time = my_toc(level_tic);
          merge_tick = merge_tick + m_time;
          merge_tick_amo = merge_tick + (m_time / chunk_merges);
          
          fprintf("\n  ** Finished: Merges for tree level %d\n ", i);
          
          % swap the variables
          U_t = U_m(1:chunk_merges, :, :);
          G_t = G_m(1:chunk_merges, :, :);
          

        end
        %fprintf("\n ** Finished Tree Merges\n");

        % final values for the subspace
        U_final = squeeze(U_t);
        % and the singular values
        G_final = squeeze(G_t);
        
        % print a nice message
        fprintf("\n -- Finised running federated trial for alpha %d\n",  ...
          p_alpha);
        
        % individual ticks
        merge_time(t, tree_depth, c_alpha) = merge_tick;
        pca_time(t, tree_depth, c_alpha) = pca_tick;
        
        % individual ticks (but amortised)
        merge_time_ta(t, tree_depth, c_alpha) = merge_tick_amo;
        pca_time_ta(t, tree_depth, c_alpha) = pca_tick / chunks;
        
        % total & amortised ticks
        t_total(t, tree_depth, c_alpha) = (pca_tick + merge_tick);
        ta_total(t, tree_depth, c_alpha) = (pca_tick + merge_tick) / chunks;
      end
      fprintf("\n -- Finished running trials for Tree Depth %d", tree_depth);
  end
  fprintf("\n -- Finished running for T=%d\n", T);
end

%%  Total Time figures (total, amortised)

% plot total worker scaling
fig = figure;
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

% print output
st = sprintf("federated_total_feat_%d_rank_%d_depth_%d", ...
  feats, rank, max_tree_depth);
print_fig(fig, st, params);

% enlarge the fonts
set(gca, 'FontSize', 18);

% plot amortised worker scaling
fig = figure;
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

% print output
st = sprintf("federated_amortised_total_feat_%d_rank_%d_depth_%d", ...
  feats, rank, max_tree_depth);
print_fig(fig, st, params);

% enlarge the fonts
set(gca, 'FontSize', 18);


%%  PCA Time figure (total)

% plot PCA time (total)
fig = figure;
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
  c_total = squeeze(pca_time(i, :, :));
  errorbar(mean(c_total, 2), std(c_total, 0, 2), '-+', 'LineWidth', 2);
  wlegs{i} = sprintf("T=%dK", (T*i) / 1000);
end
% title and stuff
title('PCA time for problem size');
ylabel('time (s)');
xlabel('node count');
xticks(1:max_tree_depth);
xticklabels(plegs);
legend(wlegs);

% print output
st = sprintf("federated_pca_time_feat_%d_rank_%d_depth_%d", ...
  feats, rank, max_tree_depth);
print_fig(fig, st, params);

% enlarge the fonts
set(gca, 'FontSize', 18);

% plot PCA time (amortised)
fig = figure;
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
  c_total = squeeze(pca_time_ta(i, :, :));
  errorbar(mean(c_total, 2), std(c_total, 0, 2), '-+', 'LineWidth', 2);
  wlegs{i} = sprintf("T=%dK", (T*i) / 1000);
end
% title and stuff
title('Amortised PCA time for problem size');
ylabel('time (s)');
xlabel('node count');
xticks(1:max_tree_depth);
xticklabels(plegs);
legend(wlegs);

% print output
st = sprintf("federated_pca_amortised_time_feat_%d_rank_%d_depth_%d", ...
  feats, rank, max_tree_depth);
print_fig(fig, st, params);

% enlarge the fonts
set(gca, 'FontSize', 18);

%% Merge time figure (total, amortised)

% plot merge ticks (total)
fig = figure;
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
  c_total = squeeze(merge_time(i, :, :));
  errorbar(mean(c_total, 2), std(c_total, 0, 2), '-+', 'LineWidth', 2);
  wlegs{i} = sprintf("T=%dK", (T*i) / 1000);
end
% title and stuff
title('Merge time per level');
ylabel('time (s)');
xlabel('node count');
xticks(1:max_tree_depth);
xticklabels(plegs);
legend(wlegs);

% print output
st = sprintf("federated_merge_time_feat_%d_rank_%d_depth_%d", ...
  feats, rank, max_tree_depth);
print_fig(fig, st, params);

% enlarge the fonts
set(gca, 'FontSize', 18);

% plot merge ticks (amortised)
fig = figure;
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
  c_total = squeeze(merge_time_ta(i, :, :));
  errorbar(mean(c_total, 2), std(c_total, 0, 2), '-+', 'LineWidth', 2);
  wlegs{i} = sprintf("T=%dK", (T*i) / 1000);
end
% title and stuff
title('Amortised Merge time per level');
ylabel('time (s)');
xlabel('node count');
xticks(1:max_tree_depth);
xticklabels(plegs);
legend(wlegs);

% print output
st = sprintf("federated_merge_amortised_time_feat_%d_rank_%d_depth_%d", ...
  feats, rank, max_tree_depth);
print_fig(fig, st, params);

% enlarge the fonts
set(gca, 'FontSize', 18);

%% Paper figure - total incl. merge and pca time)

fig = figure;

subplot(1, 3, 1);
% plot total worker scaling
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

subplot(1, 3, 2);

% plot PCA time (total)
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
  c_total = squeeze(pca_time(i, :, :));
  errorbar(mean(c_total, 2), std(c_total, 0, 2), '-+', 'LineWidth', 2);
  wlegs{i} = sprintf("T=%dK", (T*i) / 1000);
end
% title and stuff
title('PCA time for problem size');
ylabel('time (s)');
xlabel('node count');
xticks(1:max_tree_depth);
xticklabels(plegs);
legend(wlegs);


subplot(1, 3, 3);

% plot merge ticks (total)
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
  c_total = squeeze(merge_time(i, :, :));
  errorbar(mean(c_total, 2), std(c_total, 0, 2), '-+', 'LineWidth', 2);
  wlegs{i} = sprintf("T=%dK", (T*i) / 1000);
end
% title and stuff
title('Merge time per level');
ylabel('time (s)');
xlabel('node count');
xticks(1:max_tree_depth);
xticklabels(plegs);
legend(wlegs);

% print output
st = sprintf("federated_combined_total_time_feat_%d_rank_%d_depth_%d", ...
  feats, rank, max_tree_depth);
print_fig(fig, st, params);

%% Paper figure - amortised times (merge, pca)

fig = figure;


subplot(1, 3, 1);
% plot amortised worker scaling
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


subplot(1, 3, 2);

% plot PCA time (amortised)
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
  c_total = squeeze(pca_time_ta(i, :, :));
  errorbar(mean(c_total, 2), std(c_total, 0, 2), '-+', 'LineWidth', 2);
  wlegs{i} = sprintf("T=%dK", (T*i) / 1000);
end
% title and stuff
title('Amortised PCA time for problem size');
ylabel('time (s)');
xlabel('node count');
xticks(1:max_tree_depth);
xticklabels(plegs);
legend(wlegs);


subplot(1, 3, 3);

% plot merge ticks (amortised)
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
  c_total = squeeze(merge_time_ta(i, :, :));
  errorbar(mean(c_total, 2), std(c_total, 0, 2), '-+', 'LineWidth', 2);
  wlegs{i} = sprintf("T=%dK", (T*i) / 1000);
end
% title and stuff
title('Amortised Merge time per level');
ylabel('time (s)');
xlabel('node count');
xticks(1:max_tree_depth);
xticklabels(plegs);
legend(wlegs);

% print output
st = sprintf("federated_combined_amortised_time_feat_%d_rank_%d_depth_%d", ...
  feats, rank, max_tree_depth);
print_fig(fig, st, params);

%% Clean up

g_time = toc(glob_tic);

% print the total execution time
fprintf("\n ** Total Elapsed time was %d seconds\n", g_time);

% delete the parallel pool
if ~isempty(p) && destroy_pool == 1
  delete(p)
end
