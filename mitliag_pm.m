function [U, opt_out] = mitliag_pm(Y, r, params)
%MITLIAG_PM This is an implementation of the Power Method (PM). 
% As described in "Memory Limited, Streaming PCA" by Mitliagkas et al.
%
% Link: https://arxiv.org/pdf/1307.0032.pdf
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 03/06/2020
% 
% License: GPLv3
%
  fprintf('\n ** Running Mitliagkas Power Method...\n')
  
  % ambient dim defined by our input matrix Y
  [n, Ti] = size(Y);
  
  % set project data to be equal to zero
  opt_out.Yr = 0;
    
  % check if n < r or n == 1
  if n == 1 || n < r
    error(" ** ERR: Ambient dimension must be > 1 and r < n **");
  end
  
  % if we do not have params, set params to default values.
  if nargin < 3
    params.no_err = 1;
    params.no_final_err = 1;
    params.blk_size = n;
    params.use_blk_err = 1;
  % params was supplied, check if we have all entries otherwise
  % set missing to defaults.
  else
    % check if no error flag is present
    if ~isfield(params, 'no_err')
      params.no_err = 1;
    end
    
    % check if use block error flag is present
    if ~isfield(params, 'use_blk_err')
      params.use_blk_err = 1;
    end
    
    % check if no projected data flag is present
    if ~isfield(params, 'no_final_err')
      params.no_final_err = 0;
    end
    
    % check if batch block size is present
    if ~isfield(params, 'blk_size')
      params.blk_size = n;
    end
  end

  % set the block, depending on argument
  if params.blk_size < n
    B = n;
  else
    if params.blk_size < n
      fprintf("\n !! WARN: Block size must be at least n !!\n");
      B = n;
    else
      B = params.blk_size;
    end
  end
  
  % check if Ti < B, in which case we cannot run it
  if Ti < B
    error("\n Block size must be lower than the number of columns");
  end
  
  K = floor(Ti/B);              % Number of blocks
  cnt = 1;                      % index for error block align 
  
  % preallocate based on no error run and block error
  if params.no_err == 0
    if params.use_blk_err == 1
      opt_out.T = nan(1, K);        % T steps for error log
      opt_out.ErrFro = nan(1, K);   % Fro normalised error with T
    else
      opt_out.T = 1:Ti;             % T steps for error log
      opt_out.ErrFro = nan(1, Ti);  % Fro normalised error with T
    end
  % else if errors are disabled just set them to zero
  else
    opt_out.T = 0;
    opt_out.ErrFro = 0;
  end
  
  % output the number of blocks and their size
  fprintf([' ** Total number of blocks (k): %d ', ...
    'with Block size of: %d\n'], K, B);
  
  % Random initialization
  U = orth(randn(n, r));
  
  % start timing
  ts = tic;
  for k = 1:K
    % update the previous subspace estimation
    U_old = U;
    % indices for current block
    min_t = ((k-1)*B)+1;
    max_t = k*B;
    % initialise the temporary place-holder
    temp = zeros(n,r);
    for i = 1:B
      index = (k-1)*B+i;
      temp = temp+1/B*(Y(:,index)*(Y(:,index))')*U(:, 1:r);
    end 
    [U, ~, ~] = svds(temp, r);

    % Mitliagkas Power Method Errors start

    if params.no_err == 0
      % use the block error
      if params.use_blk_err == 1      
        % generate current Yr estimation
        YrHat_c = (U*U')*Y(:, 1:k*B);
        % Frobenius norm incremental error, per block at 
        % at (k-1)B:kB normalised with current T.
        temp = sum(sum((Y(:, 1:max_t)-YrHat_c).^2, 1));
        opt_out.ErrFro(cnt) = temp/max_t;
        opt_out.T(cnt) = max_t;
        cnt = cnt + 1;
      else
        % generate current Yr estimation
        YrHat_c = (U_old*U_old')*Y(:, 1:k*B);
        % Frobenius norm incremental error, per column in the block  
        % at (k-1)B:kB normalised with current T.
        temp = cumsum(sum((Y(:, 1:max_t)-YrHat_c).^2, 1));
        opt_out.ErrFro(cnt:cnt+B-1) = temp(min_t:max_t)./(min_t:max_t);
        cnt = cnt + B;
      end
    end

    % Mitliagkas Power Method errors end
    
  end
  
  % check if we require the projected data
  if params.no_err == 0
    if params.no_final_err == 0
        opt_out.Yr = (U*U')*Y(:, 1:k*B);
    end
  end
  
  % finally, set current elapsed time
  opt_out.t = my_toc(ts);   
end