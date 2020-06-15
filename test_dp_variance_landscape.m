% Comparison of streaming vs regular mod-sulq variance landscape.
%
% Based on work of Grammenos et al.: https://arxiv.org/abs/1907.08059
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 30/05/2020
% 
% License: GPLv3
%

%% Initiliase

clc; clear; close all;

%% Test variance landscape

% range for ambient dimension
dim = 50:50:1000;

% block size
blk_sz = 5:500:10000;

% privacy parameters
e_p_vals = 0.2:0.2:1.4; % better to have a range of epsilon.
delta = .1;

% compute dimensions
dim_len = size(dim, 2);
blk_len = size(blk_sz, 2);
ep_len = size(e_p_vals, 2);

% pre-allocate the arrays for speed
mod_sulq_var = zeros(ep_len, dim_len, blk_len);
smod_sulq_var = zeros(ep_len, dim_len, blk_len);

% run the loop to compute the parameters
for k = 1:ep_len
  e_p = e_p_vals(k);
  for i = 1:dim_len
    d = dim(i);
    for j = 1:blk_len
      b = blk_sz(j);
      mod_sulq_var(k, i, j) = mod_sulq_variance(d, b, e_p, delta);
      smod_sulq_var(k, i, j) = stream_mod_sulq_variance(d, b, e_p, delta);
    end
  end
end

%% Surf plots

% create the meshgrid for the landscape
[blk_axis, ambient_axis] = meshgrid(blk_sz, dim);
landscape = blk_axis.*ambient_axis;

% plot the surface landscape.
figure;
subplot(2, 1, 1);
hold on
for k = 1:ep_len
  var_val = squeeze(mod_sulq_var(k, :, :));
  surf(blk_axis, ambient_axis, var_val, landscape, 'FaceAlpha', 0.8);
  colormap(parula(k))
end
hold off
xlabel('block size');
ylabel('ambient dimension');
zlabel('variance');
title('mod sulq variance');
legend(cellstr(num2str(e_p_vals', 'e=%.1f')))
colorbar
% adjust the view point
view(50, 22)

subplot(2, 1, 2);
hold on
for k = 1:ep_len
  var_val = squeeze(smod_sulq_var(k, :, :));
  surf(blk_axis, ambient_axis, var_val, landscape, 'FaceAlpha', 0.8);
end
hold off
xlabel('block size');
ylabel('ambient dimension');
zlabel('variance');
title('stream mod sulq variance');
legend(cellstr(num2str(e_p_vals', 'e=%.1f')))
colorbar
% adjust the viewpoint
view(50, 22)


