function [params, YrHat_c] = fpca_edge_errors(Y, B, U, U_old, min_t, max_t, params)
%FPCA_EDGE_ERRORS Compute the errors of F-PCA in one function.
%
% Based on work of Grammenos et al.: https://arxiv.org/abs/1907.08059
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 30/05/2020
% 
% License: GPLv3
%

  % check if we do use errors at all
  if params.no_err == 1
    YrHat_c = 0;
    return
  end

  % normal errors start
  
  % fetch the error counter from the parameters
  err_cnt = params.err_cnt;
  
  % check if we compute block errors
  if params.use_blk_err == 1
    % Now calculate the Block normalised errors
    YrHat_c = (U*U')*Y(:, 1:max_t);
    % Frobenius norm incremental error, per block located at
    % kb normalised with current T.
    temp = sum(sum((Y(:, 1:max_t)-YrHat_c).^2, 1));
    params.ErrFro(err_cnt) = temp/max_t;
    params.T(err_cnt) = max_t;
    err_cnt = err_cnt + 1;
  else
    % Now calculate the Block normalised errors
    YrHat_c = (U_old*U_old')*Y(:, 1:max_t);
    % frobenius norm incremental error, per column in block located
    % at (k-1)b:kb normalised with current T.
    temp = cumsum(sum((Y(:, 1:max_t)-YrHat_c).^2, 1));
    params.ErrFro(err_cnt:err_cnt+params.blk_size-1) = temp(min_t:max_t)./(min_t:max_t);
    err_cnt = err_cnt + params.blk_size;
  end
  
  % update the error counter in the parameters
  params.err_cnt = err_cnt;

  % normal errors end

  % extended errors start

  if params.no_err == 0 && params.use_ext_errs == 1
    % calculate the projections, if needed
    proj_block = U'*B;

    % reconstruction energy
    cur_recon = U*proj_block;

    % update the projections block
    params.proj(1:params.r, min_t:max_t) = proj_block;

    % reconstruction add
    params.recon(:, min_t:max_t) = cur_recon;

    % add the relative errors
    params.relerrs(min_t:max_t) = (sum((B-cur_recon).^2, 1)./sum(B.^2, 1));

    % set the rpcs
    params.rpcs(min_t:max_t) = params.r;
  end

  % extended errors end
end

