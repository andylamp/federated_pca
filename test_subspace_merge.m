% The script is designed to evaluate the performance of the subspace
% merging procedure.
%
% Based on work of Grammenos et al.: https://arxiv.org/abs/1907.08059
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 02/06/2020
% 
% License: GPLv3
%

%% Initliasation
clc; clear; close all;

% for reproducibility
rng(300);

% the type used
params.type = "merge";
% enable printing
params.pflag = 0;

params = setup_vars(params);

% put everything in one plot
one_plot = 1;

fprintf("\n -- Merging Test suite starting\n\n");

%% Test execution (Basic)

fprintf("\n >> Running basic merging tests...\n");

% configuration for basic test
feats = 500;  % number of features
alpha1 = 0.1; % alpha for first distribution
alpha2 = 1;   % alpha for second distribution
T1 = 600;     % columns for first distribution
T2 = 400;     % columns for second distribution
T3 = 300;     % columns for third distribution

r = feats;  % common rank
r1 = 10;    % rank 1 test
r2 = 10;    % rank 2 test
r3 = 5;     % rank 3 test

% use basic merge to test
% 1 - naive svd merge
% 2 - naive qr merge
% 3 - block matrix merge using econ. svd and qr
algo_type = 1; 
lambda1 = 1;
lambda2 = 1;

% generate the synthetic dataset
Y1 = rand(feats, T1);
Y2 = rand(feats, T2);
Y3 = rand(feats, T3);

% perform the svds with same rank (at first)
[U1, S1, ~] = svds(Y1, r);
[U2, S2, ~] = svds(Y2, r);
[U3, S3, ~] = svds(Y3, r);
yy = [Y1, Y3, Y2];
[Uff, Sff, ~] = svds(yy, r);

% try to merge
[Uf, Sf] = fpca_subspace_merge(U1, S1, U2, S2, lambda1, lambda2, r, algo_type);
[Uf, Sf] = fpca_subspace_merge(Uf, Sf, U3, S3, lambda1, lambda2, r, algo_type);

% diff
fprintf(" ** Using equal ranks (r: %d)\n", r);
fprintf(" ** Subspace (abs) diff: %d\n", norm(abs(Uf)-abs(Uff), 'fro'));
fprintf(" ** Subspace diff: %d\n", norm(Uf-Uff, 'fro'));
fprintf(" ** Singular value diff: %d\n", norm(Sf-Sff(1:r, 1:r), 'fro'));

% perform the svds with same rank (at first)
[U1, S1, ~] = svds(Y1, r1);
[U2, S2, ~] = svds(Y2, r2);
[Uff, Sff, ~] = svds([Y1, Y2], max(r1, r2));

% try to merge
[Uf, Sf] = fpca_subspace_merge(U1, S1, U2, S2, lambda1, lambda2, max(r1, r2), algo_type);

% diff
fprintf("\n ** Using unequal ranks (r1: %d, r2: %d)\n", r1, r2);
fprintf(" ** Subspace (abs) diff: %d\n", norm(abs(Uf)-abs(Uff), 'fro'));
fprintf(" ** Subspace diff: %d\n", norm(Uf-Uff, 'fro'));
fprintf(" ** Singular value diff: %d\n", norm(Sf-Sff, 'fro'));


fprintf("\n >> Finished basic merging tests...\n");


%% Running over variable sizes to evaluate error scaling against SVD

fprintf("\n >> Running over variable sizes...\n\n");

% number of features (ambient dimension)
feats = 800;
% number of vectors (columns)
T = [feats, 2*feats, 3*feats, 4*feats, 5*feats];
% T = [200, 400, 600, 800, 1000];
% target rank
r = 100;

% synthetic dataset parameter for Power Law
synth_params.spectrum_type = "pl";
synth_params.alpha = 1;
synth_params.lambda = .01;

% preallocation of error arrays
errf_fast_u = zeros(1, size(T, 2));
errf_fast_g = zeros(1, size(T, 2));
errf_svd_u = zeros(1, size(T, 2));
errf_svd_g = zeros(1, size(T, 2));
f_times = zeros(1, size(T, 2));
s_times = zeros(1, size(T, 2));

% run for T
for i = 1:size(T, 2)
  fprintf("\n == Running for T: %d\n", T(i));
  % define the chunk size for this particular instance
  chunkSize = T(i)/2;
  % generate the data
  Y = synthetic_data_gen(feats, T(i), synth_params);
  % perform the offline r-SVD on the full dataset
  [Uf, Gf, ~] = svds(Y, r);
  
  % use halves
  
  % first half
  [Um_1, Sm_1] = fpca_edge(Y(:, 1:chunkSize), r);
  % second half
  [Um_2, Sm_2] = fpca_edge(Y(:, chunkSize+1:end), r);

  % svd merge
  s_tic = tic;
  [Um_f_svd, Sm_f_svd] = fpca_subspace_merge(Um_1, Sm_1, Um_2, Sm_2, lambda1, lambda2, r, 1);
  s_times(i) = toc(s_tic);
  
  % fast merge
  f_tic = tic;
  [Um_f_fast, Sm_f_fast] = fpca_subspace_merge(Um_1, Sm_1, Um_2, Sm_2);
  f_times(i) = toc(f_tic);
  
  % check the errors using fast method
  errf_fast_u(i) = (1/T(i)) * immse(Um_f_fast, Uf);
  errf_fast_g(i) = (1/T(i)) * immse(Sm_f_fast, Gf);
  
  % check the errors using svd method
  errf_svd_u(i) = (1/T(i)) * immse(Um_f_svd, Uf);
  errf_svd_g(i) = (1/T(i)) * immse(Sm_f_svd, Gf);
end

fprintf("\n >> Finishd running over variable sizes...");
fprintf("\n >> Plotting results.");

my_ticks = size(T, 2);

% plot U errors
figure;
subplot(1, 3, 1)
plot(1:my_ticks, errf_fast_g, '*-', 'LineWidth', 2);
hold on;
plot(1:my_ticks, errf_svd_g, '+-', 'LineWidth', 2);
hold off;
title("Errors of fast vs svd for U");
legend("fast", "svd");
xticks(1:my_ticks);
xticklabels(num2cell(T));
xlabel("T");
ylabel("error (mse)");

% plot Singular Value errors
%figure
subplot(1, 3, 2)
plot(1:my_ticks, errf_fast_u, '*-', 'LineWidth', 2);
hold on;
plot(1:my_ticks, errf_svd_u, '+-', 'LineWidth', 2);
hold off;
title("Errors of fast vs svd for Singular Values");
legend("fast", "svd");
xticks(1:my_ticks);
xticklabels(num2cell(T));
xlabel("T");
ylabel("error (mse)");

% plot exec time
%figure
subplot(1, 3, 3)
plot(f_times, '*-', 'LineWidth', 2);
hold on;
plot(s_times, '+-', 'LineWidth', 2);
hold off;
title("Time for fast vs svd merging");
legend("fast", "svd");
xticks(1:my_ticks);
xticklabels(num2cell(T));
xlabel("T");
ylabel("Time (s)");


% finally set the fonts to be larger
set(findall(gcf,'-property','FontSize'),'FontSize',14)
% make the figure larger from the get go
set(gcf, 'Units', 'Normalized', 'Position',  [.4, .1, .3, .6])

fprintf("\n -- Merging Test suite finished\n");


