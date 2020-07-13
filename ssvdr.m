function [U, S, V] = ssvdr(Y, r, params)
%SSVDR Streaming r-truncated SVD. 
%
% This can be a direct replacement to svds albeit much faster and almost 
% equivalent in terms of performance.
%
% Based on work of Grammenos et al.: https://arxiv.org/abs/1907.08059
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 14/06/2020
% 
% License: GPLv3
%

  % grab the size
  [dim, Ti] = size(Y);
  
  % check we have extra params
  if nargin < 3
    params = struct;
  end
  
  % check if we are have enabled adaptive behaviour.
  if ~isfield(params, 'adaptive')
    params.adaptive = 0;
  end
  
  % check if we are given a block, if not use a block equal to r.
  if ~isfield(params, 'blk_size')
    params.blk_size = r;
  end

  % check if we are need to return vt.
  if ~isfield(params, 'ret_vt')
    params.ret_vt = 0;
  end
  
  % do a sanity check
  if params.blk_size < r
    error("block size cannot be less than target rank.")
  end
  
  % initialise
  U = NaN;
  S = NaN;
  V = NaN;
  
  % check if we return Vt
  if params.ret_vt == 1
    V = 0;
  end
  
  % run for the blocks
  blocks = floor(Ti/params.blk_size);
  
  % check if it fits in one block.
  if blocks < 2
    if isnan(V)
      [U, S, ~] = svds(Y, r);
    else
      [U, S, V] = svds(Y, r);
    end
    return
  end
  
  % check if we have adaptive rank
  if params.adaptive == 1
    fprintf("\n ** Adaptive rank enabled.\n");
    % adjustable rank parameters
    params.dim = dim;       % ambient dimension
    params.blocks = blocks; % the number of blocks
    
    params.hcnt = 0;        % hold-off counter
    params.holdoff = 0;     % the hold-off value
    params.inc = 0;         % total rank increments
    params.dec = 0;         % total rank decrements
    params.d_rank = r;      % starting rank, equal to the r
    params.r = r;           % the current rank
    
    params.tr_lo = 1;       % low threshold mark
    params.tr_hi = 10;      % high threshold mark
    
    params.verbose = 1;     % to print rank increases/decreases

    % singular value weight percentage vector
    params.sv_weight_percentage = zeros(1, params.dim);
    % the singular value vector
    params.sv_vector = zeros(1, params.dim);
  end
  
  % otherwise run for all the blocks minus one; since, if there is spillage
  % then the last block has to be slightly larger.
  for k = 1:blocks-1
    % fetch the current B bounds from Y
    min_t = ((k-1) * params.blk_size)+1;
    max_t = k * params.blk_size;
    % slice B block from Y subject to the indices defined by min_t/max_t.
    B = Y(:, min_t:max_t);
    % update the block
    [U, S, V] = ssvdr_block_update(B, r, U, S, V);
    % check if we are adaptive
    if params.adaptive == 1
      % update the rank, if needed
      [U, S, params] = ssvdr_rank_adjust(U, S, k, params);
      r = params.r;
    end
  end
  
  % compute last min_t
  min_t = (k * params.blk_size) + 1;
  % might be one block and something, due to spillage.
  B = Y(:, min_t:end);
  % final leftovers
  [U, S, V] = ssvdr_block_update(B, r, U, S, V);
end

