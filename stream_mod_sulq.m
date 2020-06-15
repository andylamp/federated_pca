function [U, S] = stream_mod_sulq(Y, r, e_p, delta, b, thresh)
%STREAM_MOD_SULQ Perform Streaming MOD-SuLQ.
%
% Based on work of Grammenos et al.: https://arxiv.org/abs/1907.08059
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 19/04/2020
% 
% License: GPLv3
%
  
  % grab the size
  [dim, ~] = size(Y);
  
  % check if we are given a block, if not use a block equal to r.
  if nargin < 5
    b = r;
  elseif b < r
    error("block size cannot be less than target rank.")
  end
  
  % set the threshold for processing in one go (for speed).
  if nargin < 6
    thresh = 50;
  end
  
  % initialise
  U = NaN;
  S = NaN;
  
  % run for the blocks
  blocks = floor(dim/b);
  
  % check if it fits in one block.
  if blocks < 2 || dim < thresh
    % generate the mod-sulq mask for Y
    B = Y*Y';
    mask = stream_mod_sulq_mask(B, e_p, delta);
    % B = mod_sulq(Y, .1, .1);
    % now get the pc's
    [U, S, ~] = svds(B + mask, r);
    % [U, S, ~] = svds(B, r);
    return
  end
  
  % otherwise run for all the blocks minus one; since, if there is spillage
  % then the last block has to be slightly larger.
  for k = 1:blocks-1
    % fetch the current B bounds from Y
    min_t = ((k-1)*b)+1;
    max_t = k*b;
    % slice B block from Y*Y' subject to the indices defined by min_t/max_t.
    B = (1/b) * Y * Y(min_t:max_t, :)';
    % generate the mod-sulq mask for that block
    mask = stream_mod_sulq_mask(B, e_p, delta);
    % update the block
    [U, S, ~] = ssvdr_block_update(B + mask, r, U, S);
  end
  
  % compute last min_t
  min_t = (k*b) + 1;
  % might be one block and something, due to spillage.
  B = Y*Y(min_t:end, :)';
  % generate the mod-sulq mask for that block
  mask = stream_mod_sulq_mask(B, e_p, delta);
  % final leftovers
  [U, S, ~] = ssvdr_block_update(B + mask, r, U, S);
end

