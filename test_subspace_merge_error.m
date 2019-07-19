%% Test Merge Subspaces Errors
%
% Description: 
%   This code is a toy bench for testing the merge subspaces algorithm
%   and its performance along different parameters and settings.
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 18/07/2019
% 
% License: GPLv3
%

%% Initialisation
clc; clear; close all;


%% Parameter setting
feats = 400;
T = [2*feats, 3*feats, 4*feats, 5*feats];%[200, 400, 600, 800, 1000];
r = 50;

% synthetic dataset parameter for Power Law
lambda = 1;
alpha = 0.01;

% preallocation of error arrays
errf_fast_u = zeros(1, size(T, 2));
errf_fast_g = zeros(1, size(T, 2));
errf_svd_u = zeros(1, size(T, 2));
errf_svd_g = zeros(1, size(T, 2));
f_times = zeros(1, size(T, 2));
s_times = zeros(1, size(T, 2));

%% Error scaling against SVD

% run for T
for i = 1:size(T, 2)
  fprintf("\n -- Running for T: %d\n", T(i));
  % define the chunk size for this particular instance
  chunkSize = T(i)/2;
  % generate the data
  Y = synthetic_data_gen(feats, T(i), lambda, alpha);
  % perform the offline r-SVD on the full dataset
  [Uf, Gf, ~] = svds(Y, r);
  
  % use halves
  
  % SPCA
  [~, ~, Um_1, Gm_1, ~, ~, ~] = ...
    spca_edge(Y(:, 1:chunkSize), r);
  [~, ~, Um_2, Gm_2, ~, ~, ~] = ...
    spca_edge(Y(:, chunkSize:end), r);
  
  % svd merge
  s_tic = tic;
  [Um_f_svd, Gm_f_svd] = merge_subspaces(Um_1, Gm_1, Um_2, Gm_2, 1, 1, r, 3);
  s_times(i) = toc(s_tic);
  
  % fast merge
  f_tic = tic;
  [Um_f_fast, Gm_f_fast] = merge_subspaces(Um_1, Gm_1, Um_2, Gm_2);
  f_times(i) = toc(f_tic);
  
  
  % check the errors using fast method
  errf_fast_u(i) = mse(Um_f_fast, Uf);
  errf_fast_g(i) = mse(Gm_f_fast, Gf);
  
  % check the errors using svd method
  errf_svd_u(i) = mse(Um_f_svd, Uf);
  errf_svd_g(i) = mse(Gm_f_svd, Gf);
end

my_ticks = size(T, 2);

% plot U errors
figure;
plot(1:my_ticks, errf_fast_g, 'LineWidth', 2);
hold on;
plot(1:my_ticks, errf_svd_g, 'LineWidth', 2);
hold off;
title("Errors of fast vs svd for U");
legend("fast", "svd");
xticks(1:my_ticks);
xticklabels(num2cell(T));
xlabel("T");
ylabel("error (mse)");

% plot Singular Value errors
figure;
plot(1:my_ticks, errf_fast_u, 'LineWidth', 2);
hold on;
plot(1:my_ticks, errf_svd_u, 'LineWidth', 2);
hold off;
title("Errors of fast vs svd for Singular Values");
legend("fast", "svd");
xticks(1:my_ticks);
xticklabels(num2cell(T));
xlabel("T");
ylabel("error (mse)");

% plot exec time
figure;
plot(f_times, 'LineWidth', 2);
hold on;
plot(s_times, 'LineWidth', 2);
hold off;
title("Time for fast vs svd merging");
legend("fast", "svd");
xticks(1:my_ticks);
xticklabels(num2cell(T));
xlabel("T");
ylabel("Time (s)");

