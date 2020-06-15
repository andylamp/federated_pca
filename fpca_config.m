function [params] = fpca_config(dims, r, params)
%FPCA_CONFIG Function that ensures all of f-pca parameters are valid.
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 18/04/2020
% 
% License: GPLv3
%  
  
  fprintf("\n -- Configuring F-PCA");

  % set the dimensions
  params.dim = dims(1);
  params.Ti = dims(2);
  
  % check if we have an adaptive flag
  if ~isfield(params, 'adaptive')
    params.adaptive = 0;
  end
  
  % check if we have a private flag
  if ~isfield(params, 'private')
    params.private = 0; 
  end
  
  % check if we have a low threshold value
  if ~isfield(params, 'tr_lo')
    params.tr_lo = 1;
  end
  
  % check if we have a high threshold value
  if ~isfield(params, 'tr_hi')
    params.tr_hi = 10;
  end
  
  % check if we have a hold-off value
  if ~isfield(params, 'holdoff')
    params.holdoff = 0;
  end
  
  
  % check if we have a valid block size
  if ~isfield(params, 'blk_size')
    params.blk_size = r;
  end
  
  % check if we have a no error flag
  if ~isfield(params, 'no_err')
    params.no_err = 1;
  end
  
  % check if we have a extended error flag
  if ~isfield(params, 'use_ext_errs')
    params.use_ext_errs = 0;
  end
  
  % check if we have a block error flag
  if ~isfield(params, 'use_blk_err')
    params.use_blk_err = 0;
  end
  
  % check if we return the final error
  if ~isfield(params, 'no_final_err')
    params.no_final_err = 0;
  end
  
  
  % check if we have a verbose flag
  if ~isfield(params, 'verbose')
    params.verbose = 1;
  end
  
  
  % check if we have an diff. - priv. epsilon value
  if ~isfield(params, 'e_p')
    params.e_p = .4;
  end
  
  
  % check if we have an diff. - priv. delta value
  if ~isfield(params, 'delta')
    params.delta = 0.05;
  end
  
  
  % initialise hold-off counter
  params.hcnt = 0;
  
  % sets the rmax to be equal to target rank, initially.
  params.rmax = r;
  
  % set the rank
  params.r = r;

  % counter for block error align
  params.err_cnt = 1;
  
  % check if we are using verbose output
  if params.verbose == 1
    fprintf("\n\t >> Using verbose output");
  end
  
  % check if we are using adaptive rank estimation and report it.
  if params.adaptive == 1
    fprintf("\n\t >> Using Adaptive Rank Estimation");

    % adjustable rank parameters
    params.hcnt = 0;     % hold-off counter
    params.inc = 0;      % total rank increments
    params.dec = 0;      % total rank decrements
    params.d_rank = r;   % starting rank, equal to the r

    % singular value weight percentage vector
    params.sv_weight_percentage = zeros(1, params.dim);
  end                 
  
  % singular value vector
  params.sv_vector = zeros(1, params.dim);
  
  % check if we are using a differentially private F-PCA
  if params.private == 1
    fprintf("\n\t >> Using Diff. Privacy (e_p: %2.2f, delta: %2.2f)", ...
      params.e_p, params.delta) 
%     
%     % compute the number of blocks for streaming mod-sulq
%     if params.blk_size > params.dim
%       params.priv_blocks = floor(params.dim/2);
%     else
%       params.priv_blocks = floor(params.d/params.blk_size);
%     end
  end

  % check we calculate the error (disabled for speed runs)
  if params.no_err == 0
    fprintf("\n\t >> Using error computation (slow)");
    
    % check if we have raised the block error flag
    if params.use_blk_err == 1
      fprintf("\n\t  == Using block-based error computation");
    else
      fprintf("\n\t  == Using analytical error computation");
    end
    
    % check if we have raised the extended error flag
    if params.use_ext_errs == 1
      fprintf("\n\t  == Using extended error computation");
    else
      fprintf("\n\t  == Not using extended error computation");
    end
  else
    fprintf("\n\t >> Skipping error computation (fast)");
    % errors, if needed
    params.ErrFro = 0;
    params.T = 0;

    % extensive errors (proj., rel., and others).
    params.use_ext_errs = 0;
  end

  % check if n < r or n == 1
  if params.dim == 1 || params.dim < r
    error("\n ** ERR: Ambient dimension must be > 1 and r < n **");
  elseif params.verbose == 1
    fprintf("\n\t >> Processing input of %d x %d", params.dim, params.Ti);
  end
  
  % check hold-off value
  if params.holdoff < 0
    error(" ** ERR: hold-off value of must be greater or equal to one");
  elseif params.verbose == 1
    % show hold-off ticks
    fprintf("\n\t >> hold-off value: %d", params.holdoff);
  end

  % check the cut-off thresholds
  if params.tr_lo > params.tr_hi
    error("\n ** ERR: low threshold greater than high");
  elseif params.verbose == 1
    fprintf("\n\t >> Using cut-off thresholds, low: %2.2f and high %2.2f", ...
       params.tr_lo, params.tr_hi);
  end
  
  % check the block size
  if params.blk_size < r
    fprintf(['\n !! WARN: Block size must be at least r,', ...
      ' resetting to default b=r !!\n']);
    params.blk_size = r;
  elseif params.Ti < params.blk_size
    error("\n ** ERR: Block size must be lower than the number of columns");
  end
  
  % compute number of blocks
  params.blocks = floor(params.Ti/params.blk_size);
  if params.verbose == 1  
    % output the block number
    fprintf(['\n\t >> Number of blocks to process: %d ', ...
      'using block size of: %d'], params.blocks, params.blk_size);
  end
  
  % preallocate error matrices
  %
  % T: steps for error log
  % ErrFro: Fro normalised error with T
  if params.no_err == 0
    % different matrice sizes based on analytical vs block error.
    if params.use_blk_err == 1
      params.T = nan(1, params.blocks);       
      params.ErrFro = nan(1, params.blocks);  
    else
      params.T = 1:params.Ti;  
      params.ErrFro = nan(1, size(params.T, 2));
    end
    
    % use extended error matrices, for relative errors, pc's evolution,
    % and so on - if not tracked, they are set to NaN.
    if params.use_ext_errs == 1      
      params.recon = zeros(params.dim, params.Ti);
      params.rpcs = zeros(1, params.Ti);
      params.relerrs = zeros(1, params.Ti);
      params.proj = zeros(params.dim, params.Ti);
    else
      params.recon = NaN;
      params.rpcs = NaN;
      params.relerrs = NaN;
      params.proj = NaN;      
    end
  end
  
  % knock for finish
  fprintf("\n -- Finished configuring F-PCA\n");
  
end

