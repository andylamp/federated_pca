function [U, S, params] = ssvdr_rank_adjust(U, S, k, params)
%SSVDR_RANK_ADJUST Performs the rank adjust of ssvdr.
%
% This is supplemendary to ssvdr.m.
%
% Based on work of Grammenos et al.: https://arxiv.org/abs/1907.08059
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 14/06/2020
% 
% License: GPLv3
%

  % increase the hold off counter.
  params.hcnt = params.hcnt + 1;
  
  % assign the singular value vectors
  if ~isvector(S)
    params.sv_vector(1:params.r) = diag(S(1:params.r, 1:params.r))';
  else
    params.sv_vector(1:params.r) = S(1:params.r);
  end
  
  % get the singular value sum.
  sg_sum = sum(params.sv_vector(1:params.r));
  
  % compute the percentages.
  params.sv_weight_percentage(1:params.r) = (params.sv_vector(1:params.r) ./ sg_sum) * 100;
  
  % compute the weights for first and last singular value.
  first_sg_weight = (params.sv_vector(1) / sg_sum) * 100;
  last_sg_weight = (params.sv_vector(params.r) / sg_sum) * 100;
    
  % check if we can rank adjust.
  if(params.hcnt >= params.holdoff && k < params.blocks-1)
    % check if we need to adjust - also: need to keep r < b.
    if(last_sg_weight > params.tr_hi && params.r < params.blk_size) % r < dim
      params.inc = params.inc + 1;
      % incr = 1;
      params.d_rank = params.d_rank + 1;
      params.r = params.d_rank;

      % set the r-th canonical vector.
      pcs = zeros(params.dim, 1);
      pcs(params.r+1, 1) = 1;

      % append to subspace.
      U = [U pcs];

      % add it to the lower part and make it orthonormal!
      U = qr(U, 0);
      % reset the hold-off counter.
      params.hcnt = 0;

      % adjust the S, after the increase.
      S = params.sv_vector(1:params.r);
      
      % check if we are verbose to report the change.
      if params.verbose == 1
        fprintf(['\t !! Increasing r to %d (of %d) at block %d', ...
          ' (of %d) with singular value weight: %3.2f %%', ...
          ' (first %3.2f %%)\n'], ...
          params.r, params.dim, k, params.blocks, last_sg_weight, first_sg_weight);
      end
    % check if we can adjust for lower.
    elseif (last_sg_weight < params.tr_lo && params.r > 1) % r > rank_seed
      params.d_rank = params.d_rank - 1;
      % remove a pc (no need to re-orth. since it already is).
      U = U(:, 1:(end-1));
      params.r = params.r - 1;

      % increase the decrement counter.
      params.dec = params.dec + 1;
      % reset hold-off counter.
      params.hcnt = 0;

      % adjust the S, after the increase.
      S = params.sv_vector(1:params.r);

      % check if we are verbose to report the change.
      if params.verbose == 1
        fprintf(['\t !! Decreasing r to %d (of %d) at block %d ', ...
          '(of %d) with last singular value weight: %3.2f %%', ...
          ' (first %3.2f %%)\n'], ...
          params.d_rank, params.dim, k, params.blocks, last_sg_weight, first_sg_weight);
      end
      %
    end
    %
  end
  %
end

