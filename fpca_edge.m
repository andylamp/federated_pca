function [U, S, opt_out] = fpca_edge(Y, r, params)
%FPCA_EDGE Runs the local updates for the federated FPCA scheme.
%
% Based on work of Grammenos et al.: https://arxiv.org/abs/1907.08059
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 30/05/2020
% 
% License: GPLv3
%
  fprintf('\n -- Running F-PCA Edge...');

  % configuration of fpca
  if nargin < 3
    params = fpca_config(size(Y), r);
  else
    params = fpca_config(size(Y), r, params);
  end
  
  % set the block
  b = params.blk_size;
  % compute the number of blocks
  blocks = params.blocks;              
  
  % start timing
  ts = tic;

  % initialise
  U = NaN;
  S = NaN;
  Yest = 0;
  sv_diag = zeros(1, params.dim);
  
  % run for the remaining blocks
  for k = 1:blocks
    % fetch the current y_t from Y
    min_t = ((k-1)*b)+1;
    max_t = k*b;
    B = Y(:, min_t:max_t);
    
    % check if we are doing the private version, so we handle blocks
    % differently.
    if params.private == 1
      % compute the covariance matrix in a streaming way
      [Ub, ~] = fpca_stream_mod_sulq(B, params.r, params.e_p, ...
        params.delta, b, params.Ti);
      % update the global block
      [U, S, ~] = ssvdr_block_update(Ub, params.r, U, sv_diag(1:params.r));
    else
      % run normal block update
      [U, S, ~] = ssvdr_block_update(B, params.r, U, sv_diag(1:params.r));
    end
    
    % set the singular values
    sv_diag(1:params.r) = diag(S(1:params.r, 1:params.r));
    params.sv_vector(1:params.r) = sv_diag(1:params.r);
    
    % compute errors
    [params, Yest] = fpca_edge_errors(Y, B, U, U, min_t, max_t, params);
    
    % adjust rank
    [U, params] = fpca_rank_adjust(U, k, max_t, params);
  end
  
  % just an aggregate output at the end.
  sg_sum = sum(params.sv_vector(1:params.r));
  first_sg_weight = (params.sv_vector(1) / sg_sum) * 100;
  last_sg_weight = (params.sv_vector(params.r) / sg_sum) * 100;
  fprintf(['\t !! Final last singular value weight: %3.2f %% ', ...
    '(first: %3.2f %%) using %d PCs.\n'], ...
    last_sg_weight, first_sg_weight, params.r);
  
  % check if are adaptive - if so report brief stats about it.
  if params.adaptive == 1
    fprintf("\t >> Total decs: %d, Total incs: %d, and Max rank: %d.\n", ...
      params.dec, params.inc, params.rmax)
  end
  
  % finally update finalised the r-svds estimates
  U = U(:, 1:params.r);
  S = S(1:params.r, 1:params.r);
  
  % calculate the current trial execution delta
  t = my_toc(ts);
  
  % assign these parameters only if errors are tracked
  if params.no_err == 0
    opt_out.T = params.T;
    opt_out.ErrFro = params.ErrFro;
    opt_out.recon = params.recon;
    opt_out.rmax = params.rmax;
    opt_out.rpcs = params.rpcs;
    opt_out.relerrs = params.relerrs;
  
    % also set the final error estimate of Y
    if params.no_final_err == 0
      opt_out.Yr = Yest;  
    end
  end
  
  % set the execution time
  opt_out.t = t;
  
  % knock for finish
  fprintf(' -- Finished F-PCA Edge...\n\n');
end
