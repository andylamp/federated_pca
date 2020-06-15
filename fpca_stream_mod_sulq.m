function [U, S] = fpca_stream_mod_sulq(B, r, e_p, delta, b, gTi)
%FPCA_STREAM_MOD_SULQ F-PCA Streaming MOD-SuLQ.
%
% Based on work of Grammenos et al.: https://arxiv.org/abs/1907.08059
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 30/05/2020
% 
% License: GPLv3
% 
  % grab the size
  [d, Ti] = size(B);
  
  % check if we are given a block, if not use a block equal to r.
  if nargin < 5
    b = 2*r;
  elseif b < r
    error("block size cannot be less than target rank.")
  end
  
  % initialise
  U = NaN;
  S = NaN;
  V = NaN;
  
  % run for the blocks
  blocks = floor(d/b);
  
  % check if it fits in one block.
  if blocks < 2
    % the stream MOD-SuLQ perturbation mask
    mask = stream_mod_sulq_mask([d, d], gTi, e_p, delta);
    [U, S, ~] = svds(((1/Ti) * (B*B')) + mask, r);
    return
  end
  
  % otherwise run for all the blocks minus one; since, if there is spillage
  % then the last block has to be slightly larger.
  for k = 1:blocks-1
    % fetch the current B bounds from Y
    min_t = ((k-1)*b)+1;
    max_t = k*b;
    % slice B block from Y subject to the indices defined by min_t/max_t.
    Bc = (1/Ti) * (B * B(min_t:max_t, :)');
    % the stream MOD-SuLQ perturbation mask
    mask = stream_mod_sulq_mask(size(Bc), gTi, e_p, delta);
    % update the block
    [U, S, V] = ssvdr_block_update(Bc + mask, r, U, S, V);
  end
  
  % compute last min_t
  min_t = (k*b) + 1;
  % slice B block from Y subject to the indices defined by min_t/max_t.
  Bc = (1/Ti) * (B * B(min_t:end, :)');
  % the stream MOD-SuLQ perturbation mask
  mask = stream_mod_sulq_mask(size(Bc), gTi, e_p, delta);
  % update the block
  [U, S, ~] = ssvdr_block_update(Bc + mask, r, U, S, V);
end

