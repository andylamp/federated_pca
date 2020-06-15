function [Vk, A] = mod_sulq(block, e_p, delta, symmetric, scaler)
%MOD_SULQ An implementation of Mod-SuLQ.
%
% Based on work of Chaudhuri et al.: https://arxiv.org/pdf/1207.2812
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 30/05/2020
% 
% License: GPLv3
%
  % get the block
  [d, n] = size(block);
  
  % mean
  mu = 0;
  
  % check if we apply symmetric mod-sulq
  if nargin < 4
    symmetric = 1;
  end
  
  % check if we have a scaler for the variance (wrt to number of columns)
  if nargin < 5
    scaler = 1;
  end
  
  % get the variance based on the given parameters
  beta = mod_sulq_variance(d, scaler*n, e_p, delta);
  
  % expand the covariance matrix
  A = (1/n) * (block*block');
  
  % generate the the perturbation mask
  N = mu + beta * randn(d);
  % convert N into a symmetric matrix
  if symmetric == 1
    N = N - tril(N, -1) + triu(A, 1)';
  end
  
  % final matrix
  Vk = A + N;
  
end