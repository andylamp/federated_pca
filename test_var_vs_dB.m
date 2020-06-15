% This script is designed to test the variability of the variance based
% on different rations of the ambient dimension (d) over the block size
% (B). This can provide useful insights on how this evolves as d approaches
% the dataset size.
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

clear; close all; clc;

% columns in the dataset
Ti = [1, 100:100:1000];
% the ambient dimension / data size ratio (d/B)
dB_ratio = 0:.001:0.1;
dB_ratio = dB_ratio(2:end);

% privacy parameters
e_p = 0.4;
delta = .1;

% preallocate arrays for speed
db_len = size(dB_ratio, 2);
cols_len = size(Ti, 2);
beta_xn_t = zeros(db_len, cols_len);

%% Test for different d/B ratios

% Run the loop to compute the parameters
for i = 1:cols_len
  n = Ti(i);
  for j = 1:db_len
    x = dB_ratio(j);
    beta_xn_t(i, j) = beta_xn(x, n, e_p, delta); 
  end
end

%% Plot the curves

% plot the db curve
figure;
hold on
for k = 1:db_len
  beta_val = squeeze(beta_xn_t(k, :));
  % xn_val = squeeze(cols(k, :));
  plot(beta_val);
end
legend(cellstr(num2str(Ti', 'n=%.1f')), 'location', 'best');

% get xticks len
xticks_vals = xticks;
xticks_len = size(xticks_vals, 2);

xticks_labels = cell(1, xticks_len);
xticks_labels{1} = xticks_vals(1);

for i = 2:xticks_len
  xticks_labels{i} = dB_ratio(xticks_vals(i));
end
% finally set them up
xticklabels(xticks_labels);

title('variance over d/B ratio')
xlabel('d/B ratio');
ylabel('variance value');
