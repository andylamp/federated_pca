function [T, ErrFro, Sk, Gk, Qk, Yr, t] = spca_edge(Y, r, blk_size, ...
  floor_mul, no_err) 
%SPCA_EDGE function is responsible for computing SPCA on "edge" nodes.
%
% Notes:
%  - It has to be noted that this is version of the algorithm that has the 
%    ability to calculate a lot of error metrics and flags used for its 
%    evaluation. Hence, should not be used in production as both 
%    computational and memory costs are increased.
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 18/07/2019
% 
% License: GPLv3
%
  % scope in global variables
  global use_blk_err
  global allow_print
  
  if allow_print == 1
    fprintf('\n ** Running SPCA Edge...\n');
  end
  
  % get Y details
  [dim, Ti] = size(Y);            
  
  % check we calculate the error (disabled for speed runs)
  if nargin < 5
    no_err = 0;
  end
  
  % check if n < r or n == 1
  if dim == 1 || dim < r
    error(" ** ERR: Ambient dimension must be > 1 and r < n **");
  end

  % check if we have a floor multiplier as an argument
  if nargin < 4
    floor_mul = 2;
  end

  % moses configuration
  
  % set the block, depending on argument
  if nargin < 3
    b = r;
  else
    if blk_size < r
      fprintf(['\n !! WARN: Block size must be at least r,', ...
        ' resetting to default b=r !!\n']);
      b = r;
    else
      b = blk_size;
    end
  end
  
  % check if Ti < b, in which case we cannot run it
  if Ti < b
    error("\n Block size must be lower than the number of columns");
  end

  K = floor(Ti/b);              % Number of blocks          
  cnt = 1;                      % counter for block error align
  
  % preallocate based on no error run and block error
  if no_err == 0
    if use_blk_err == 1
      T = nan(1, K);                % T steps for error log
      ErrFro = nan(1, K);           % Fro normalised error with T
    else
      T = 1:Ti;                     % T steps for error log
      ErrFro = nan(1, size(T, 2));  % Fro normalised error with T
    end
  else
    T = 0;
    ErrFro = 0;
  end

  % output the block number
  if allow_print == 1
    fprintf([' ** Total number of blocks (k): %d ', ...
      'with Block size of: %d\n'], K, b);
  end
  
  % start timing
  ts = tic;
  % form y_k, which comprises the first block
  y_k = Y(:, 1:b);
  % get the first estimation of the r-svds
  [S_k, G_k, q_k] = svds(y_k, floor(floor_mul * r));
  % reduce the dimension of the svds components
  S_k = S_k(:, 1:r);
  G_k = G_k(1:r, 1:r);
  q_k = q_k(:, 1:r);
  
  % preallocate the projections
  Proj = zeros(r, Ti);
  
  % calculate the first block projections
  proj_block = S_k'*y_k;
  
  % energy initialisation
  
  % reconstruction energy
  enSumYr = sum(sum((S_k*proj_block).^2));
  % actual energy
  enSumY = sum(sum(y_k.^2));
  
  % assign projection block
  Proj(:, 1:b) = proj_block;
  
  % errors start

  if no_err == 0
    if use_blk_err == 1
      % Now calculate the Block normalised errors
      YrHat_c = (S_k*S_k')*Y(:, 1:b);
      % Frobenius norm incremental error, per block located at
      % kb normalised with current T.
      temp = sum(sum((Y(:, 1:b)-YrHat_c).^2, 1));
      ErrFro(cnt) = temp/b;
      T(cnt) = b;
      cnt = cnt + 1;
    else
      % Now calculate the Block normalised errors
      YrHat_c = (S_k*S_k')*Y(:, 1:b);
      % frobenius norm incremental error, per column in block located
      % at (k-1)b:kb normalised with current T.
      temp = cumsum(sum((Y(:, 1:b)-YrHat_c).^2, 1));
      ErrFro(cnt:cnt+b-1) = temp(1:b)./(1:b);
      cnt = cnt + b;
    end
  end

  % errors end
  
  % run for the remaining blocks
  for k = 2:K
    S_kold = S_k;
    % fetch the current y_t
    min_t = ((k-1)*b)+1;
    max_t = k*b;
    y_k = Y(:, min_t:max_t);
    % construct the q_k
    q_k = S_k'*y_k;
    % construct the z_k
    z_k = y_k - S_k*q_k;
    % get the (economy) QR of z_k
    [s_k, v_k] = qr(z_k, 0);
    % now construct the following block matrix as is shown in the
    % algorithm in our paper:
    %
    %           |      G_k      q_k |
    % blk_mat = |                   |
    %           | zeros(zr, r)  v_k |
    %
    zr = min(b, size(v_k, 1));
    %cc_top = [G_k(1:r, :), q_k]; size(cc_top)
    %cc_bot = [zeros(zr, r), v_k]; size(cc_bot)
    %return
    blk_mat_k = [ G_k, q_k; zeros(zr, r), v_k ];
    
    % now take the r-svds of that matrix
    [u_k, G_k, q_k] = svds(blk_mat_k, floor(floor_mul * r));
    % reduce the dimension of the svds components
    u_k = u_k(:, 1:r);
    G_k = G_k(1:r, 1:r);
    q_k = q_k(:, 1:r);
    % now update the actual S_k estimation
    S_k = [S_k, s_k]*u_k;
    
    % errors start

    if no_err == 0
      if use_blk_err == 1
        % Now calculate the Block normalised errors
        YrHat_c = (S_k*S_k')*Y(:, 1:max_t);
        % Frobenius norm incremental error, per block located at 
        % kb normalised with current T.
        temp = sum(sum((Y(:, 1:max_t)-YrHat_c).^2, 1));
        ErrFro(cnt) = temp/max_t;
        T(cnt) = max_t;
        cnt = cnt + 1;
      else
        % Now calculate the Block normalised errors
        YrHat_c = (S_kold*S_kold')*Y(:, 1:max_t);
        % frobenius norm incremental error, per column in block located  
        % at (k-1)b:kb normalised with current T.
        temp = cumsum(sum((Y(:, 1:max_t)-YrHat_c).^2, 1));
        ErrFro(cnt:cnt+b-1) = temp(min_t:max_t)./(min_t:max_t);
        cnt = cnt + b;
      end
    end

    % errors end
    
    % Adjust rank start
    
    % calculate the projections, if needed
    proj_block = S_k'*y_k;
    
    % update reconstruction energy
    enSumYr = enSumYr + sum(sum((S_k*proj_block).^2));
    % update the actual energy
    enSumY = enSumY + sum(sum(y_k.^2));
    
    % update the projections block
    Proj(:, min_t:max_t) = proj_block;
    % calculate the energy ratio
    %enRatio = 100*enSumYr/enSumY;
    
    
  end
  % finally update finalised the r-svds estimates
  Sk = S_k;
  Gk = G_k;
  Qk = q_k;
  % also set the final estimate of Y
  if no_err == 0
    Yr = YrHat_c;     % update the final Yr estimation
  else
    Yr = (S_k*S_k')*Y(:, 1:max_t);
  end
  % calcualte the current trial execution delta
  t = my_toc(ts, allow_print);
  % knock for finish
  if allow_print == 1
    fprintf(' ** Finished running SPCA Edge...\n\n');
  end
end

