function [T, ErrFro, Sk, Gk, Qk, Yr, t, recon, rmax, relerrs, rpcs] = ...
  sapca_edge(Y, rank_seed, tr_bounds, holdoff, blk_size, floor_mul, ...
  no_err, silent)
%SAPCA_EDGE function is responsible for computing SAPCA on "edge" nodes.
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

  fprintf('\n ** Running SAPCA Edge...\n');

  % scope in global variables
  global use_blk_err

  % get Y details
  [dim, Ti] = size(Y);

  % get the r seed
  r = rank_seed;

  recon = zeros(dim, Ti);
  rpcs = zeros(1, Ti);
  relerrs = zeros(1, Ti);

  if nargin < 8
    silent = 1;
  end

  % check we calculate the error (disabled for speed runs)
  if nargin < 7
    no_err = 0;
  end

  % check if n < r or n == 1
  if dim == 1 || dim < r
    error(" ** ERR: Ambient dimension must be > 1 and r < n **");
  end

  % check if we have a floor multiplier as an argument
  if nargin < 6
    floor_mul = 2;
  end

  % configuration

  % set the block, depending on argument
  if nargin < 5
    b = r;
    %b = 2*r
  else
    if blk_size < r
      fprintf(['\n !! WARN: Block size must be at least r,', ...
        ' resetting to default b=r !!\n']);
      b = r;
    else
      b = blk_size;
    end
  end

  % set the hold-off value for change
  if nargin < 4
    holdoff = 0;  % hold-off blocks
  end

  % set the weight percentage threshold for change
  if nargin < 3
      hi_tr = 10; % threshold for upper weight percentage bound
      lo_tr = 1;  % threshold for lower weight percentage bound
  else
      hi_tr = tr_bounds(2);
      lo_tr = tr_bounds(1);
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
  fprintf([' ** Total number of blocks (k): %d ', ...
    'with Block size of: %d\n'], K, b);

  % adjustable rank parameters
  rmax = r;     % current max rank
  hcnt = 0;     % hold-off counter
  inc = 0;      % total rank increments
  dec = 0;      % total rank decrements
  d_rank = r;   % starting rank, equal to the block
  %incr = min(20, b);     %min(b, 20);     % increment rate
  %incr_hi = min(b, 1);
  %incr_lo = min(b, 1);
  % lambdas for bias
  gw_dd = zeros(1, dim);
  sg_weight_percentage = zeros(1, dim);

  % start timing
  ts = tic;
  % form y_k, which comprises the first block
  y_k = Y(:, 1:b);
  % get the first estimation of the r-svds
  [S_k, G_k, q_k] = svds(y_k, floor(floor_mul * r));
  % set the singular values
  gk_diag = diag(G_k);
  gw_dd(1:r) = gk_diag(1:r);

  % reduce the dimension of the svds components
  %S_k = S_k_f(:, 1:r);
  %G_k = G_k_f(1:r, 1:r);
  %q_k = q_k_f(:, 1:r);

  % preallocate the projections
  Proj = zeros(dim, Ti);

  % calculate the first block projections
  proj_block = S_k(:, 1:r)'*y_k;

  % energy initialisation

  % reconstruction energy
  cur_recon = S_k(:, 1:r)*proj_block;
  %enSumYr = sum(sum(cur_recon.^2, 1));
  % actual energy
  %enSumY = sum(sum(y_k.^2, 1));

  relerrs(1:b) = sum((y_k - cur_recon).^2, 1)./sum(y_k.^2, 1);

  % assign projection block
  Proj(1:r, 1:b) = proj_block;

  % assign the reconstruction
  recon(:, 1:b) = cur_recon;

  % set the rpcs in the first block
  rpcs(1:b) = r;
  
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

    S_kold = S_k(:, 1:r);
    % fetch the current y_t
    %k
    min_t = ((k-1)*b)+1;
    max_t = k*b;
    y_k = Y(:, min_t:max_t);

    % construct the q_k (projections)
    q_k = S_k(:, 1:r)'*y_k;
    % reconstruction energy
    cur_recon = S_k(:, 1:r)*q_k;
    % construct the z_k
    z_k = y_k - cur_recon;

    % construct the q_k (projections)
    %q_k = S_k(:, 1:r)'*y_k;
    % construct the z_k
    %z_k = y_k - S_k(:, 1:r)*q_k;
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
    % grab the previous diagonal eigenvalues
    G_k = diag(gw_dd(1:r));
    %cc_top = [G_k(1:r, 1:r), q_k(:, :)]; %size(cc_top)
    %cc_bot = [zeros(zr, r), v_k]; %size(cc_bot)
    blk_mat_k = [ G_k(1:r, 1:r), q_k(:, :); zeros(zr, r), v_k ];

    % now take the r-svds of that matrix
    %[u_k, G_k, q_k] = svds(blk_mat_k, floor(floor_mul * r));
    [u_k, G_k, q_k] = svds(blk_mat_k, floor(floor_mul * r));
    % set the singular values
    gk_diag = diag(G_k);
    gw_dd(1:r) = gk_diag(1:r);
    % reduce the dimension of the svds components
    %u_k = u_k(:, 1:r);
    %G_k = G_k(1:r, 1:r);
    %q_k = q_k(:, 1:r);
    % now update the actual S_k estimation
    %tS_k = S_k;
    blk_st = [S_k(:, 1:r), s_k];
    u_k_sz = size(u_k, 1);
    S_k = blk_st*u_k(1:u_k_sz, 1:r);

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

    % reconstruction energy
    cur_recon = S_k*proj_block;
    % update reconstruction energy
    %enSumYr = lam_fc*enSumYr + sum(sum((cur_recon).^2));
    % update the actual energy
    %enSumY = lam_fc*enSumY + sum(sum(y_k.^2, 1));

    % update the projections block
    Proj(1:r, min_t:max_t) = proj_block;

    % reconstruction add
    recon(:, min_t:max_t) = cur_recon;

    % add the relative errors
    relerrs(min_t:max_t) = (sum(z_k.^2, 1)./sum(y_k.^2, 1));

    % calculate the energy ratio
    %enRatio = 100*enSumYr/enSumY;

    % increase hold off counter
    hcnt = hcnt + 1;

    % the solution is hiding with eigenvalues!
    sg_sum = sum(gw_dd(1:r));
    sg_weight_percentage(1:r) = (gw_dd(1:r) ./ sg_sum) * 100;
    first_sg_weight = (gw_dd(1) / sg_sum) * 100;
    last_sg_weight = (gw_dd(r) / sg_sum) * 100;

    adaptable = 1;

    % set the rpcs in the block
    rpcs(min_t:max_t) = r;
    % good (with last !! high: 4, low: 2

    % check if we need to adjust
    if(adaptable == 1 && last_sg_weight > hi_tr && r < dim && hcnt >= holdoff) % r < dim
      inc = inc+1;
      incr = 1;
      d_rank = d_rank + incr;
      r = d_rank;
      pcs = zeros(dim, incr);

      for i = 1:incr
        pcs(r+i, i) = 1;
      end
      %size(S_k)
      %size(pcs)
      S_k = [S_k pcs];
      %S_k = [tS_k(:, 1:r), s_k]*u_k;
      % add it to the lower part and make it orthonormal!
      S_k = qr(S_k, 0);
      hcnt = 0;

      if silent == 0
        fprintf(['\t !! Increasing r to %d (of %d) at block %d', ...
          ' (of %d) with singular value weight: %6.2f %%', ...
          ' (first %6.2f %%)\n'], ...
          r, dim, k, K, last_sg_weight, first_sg_weight);
      end
    %  enSumYr > hi_tr*enSumY
    elseif (adaptable == 1 && last_sg_weight < lo_tr && r > 1 && hcnt >= holdoff) % r > rank_seed
      d_rank = d_rank-1;
      % remove a pc and re-orthonormalise
      %S_k = qr(S_k(:, 1:(end-1)), 0);
      S_k = S_k(:, 1:(end-1));
      r = r - 1;
      dec = dec+1;
      hcnt = 0;
      if silent == 0
        fprintf(['\t !! Decreasing r to %d (of %d) at block %d ', ...
          '(of %d)[with last singular value weight: %6.2f %%', ...
          ' (first %6.2f %%)\n'], ...
          d_rank, dim, k, K, last_sg_weight, first_sg_weight);
      end

    else
      %fprintf('\t -- Energy ratio at block %d (of %d) with r %d is: %6.2f %%\n', ...
      %  k, K, d_rank, enRatio);
    end

    % set the pcs based on new value
    rpcs(max_t) = r;
    % check for max rank
    if(rmax < r)
      rmax = r;
    end

    % Adjust rank end

  end
  fprintf(['\t !! Final last singular value weight: %6.2f %% ', ...
    '(first: %6.2f %%) using %d PCs [total decs: %d, ', ...
    'total incs: %d, max r: %d]\n'], ...
    last_sg_weight, first_sg_weight, r, dec, inc, rmax);
  % finally update finalised the r-svds estimates
  Sk = S_k(:, 1:r);
  Gk = G_k(1:r, 1:r);
  Qk = q_k(:, 1:r);
  % also set the final estimate of Y
  if no_err == 0
    Yr = YrHat_c;     % update the final Yr estimation
  else
    Yr = (S_k*S_k')*Y(:, 1:max_t);
  end
  % calculate the current trial execution delta
  t = my_toc(ts);
  % knock for finish
  fprintf(' ** Finished running SAPCA Edge...\n\n');
end
