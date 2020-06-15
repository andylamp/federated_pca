function [U, params] = fpca_rank_adjust(U, k, max_t, params)
%FPCA_RANK_ADJUST Perform rank adjust for F-PCA.
%
% Based on work of Grammenos et al.: https://arxiv.org/abs/1907.08059
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 30/05/2020
% 
% License: GPLv3
% 
    % increase hold off counter
    params.hcnt = params.hcnt + 1;
    
    % Adjust rank start
    
    % if we do not adjust
    if params.adaptive == 0
      return
    end

    % sum the singular values to get the normaliser
    sg_sum = sum(params.sv_vector(1:params.r));
    % compute the percentages
    params.sv_weight_percentage(1:params.r) = (params.sv_vector(1:params.r) ./ sg_sum) * 100;
    
    % compute the weights for first and last singular value
    first_sg_weight = (params.sv_vector(1) / sg_sum) * 100;
    last_sg_weight = (params.sv_vector(params.r) / sg_sum) * 100;
    
    % check if we can rank adjust
    if(params.hcnt >= params.holdoff && k < params.blocks-1)
      % check if we need to adjust - also: need to keep r < b.
      if(last_sg_weight > params.tr_hi && params.r < params.blk_size) % r < dim
        params.inc = params.inc + 1;
        % incr = 1;
        params.d_rank = params.d_rank + 1;
        params.r = params.d_rank;
        
        % set the r-th canonical vector
        pcs = zeros(params.dim, 1);
        pcs(params.r+1, 1) = 1;
        
        % append to subspace
        U = [U pcs];
        
        % add it to the lower part and make it orthonormal!
        U = qr(U, 0);
        % reset the hold-off counter
        params.hcnt = 0;

        % check if we are verbose to report the change.
        if params.verbose == 1
          fprintf(['\t !! Increasing r to %d (of %d) at block %d', ...
            ' (of %d) with singular value weight: %3.2f %%', ...
            ' (first %3.2f %%)\n'], ...
            params.r, params.dim, k, params.blocks, last_sg_weight, first_sg_weight);
        end
      % check if we can adjust for lower
      elseif (last_sg_weight < params.tr_lo && params.r > 1) % r > rank_seed
        params.d_rank = params.d_rank - 1;
        % remove a pc (no need to re-orth. since it already is).
        U = U(:, 1:(end-1));
        params.r = params.r - 1;
        
        % increase the decrement counter
        params.dec = params.dec + 1;
        % reset hold-off counter
        params.hcnt = 0;
        
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
    
    % set the pcs based on new value, if we track it
    if params.use_ext_errs == 1
      params.rpcs(max_t) = params.r;
    end
    
    % check for max rank
    if(params.rmax < params.r)
      params.rmax = params.r;
    end

    % Adjust rank end

end

