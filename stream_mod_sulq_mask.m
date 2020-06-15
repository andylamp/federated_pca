function [mask] = stream_mod_sulq_mask(dims, c, e_p, delta, scaler)
%STREAM_MOD_SULQ_MASK Generate the perturbation mask for block B.
%
% Function that generates the required perturbation mask per block
% subject to the correct variance (omega). There is also the option to set 
% a variance scaler if desired.
%
% Based on work of Grammenos et al.: https://arxiv.org/abs/1907.08059
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 19/04/2020
% 
% License: GPLv3
% 
  % ambient dimension (d)
  d = dims(1);
  % the block size (b)
  b = dims(2);
  
  % check if we have a specific variance scaler
  if nargin < 5
    scaler = 1;
  end
  
  % generate the omega variance
  % omega = stream_mod_sulq_variance(d, c, e_p, delta, scaler);
  omega = stream_mod_sulq_variance(d, c, e_p, delta, scaler);
  % omega = mod_sulq_variance(d, c, e_p, delta);
  
  % apply variance to random values
  % mask = (scaler * omega) * randn(d, c);
  mask = (scaler * omega) * randn(d, b);
end

