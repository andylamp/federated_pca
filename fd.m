function [Bout, opt_out] = fd(Y, ell, params)
%FD An implementation of the Frequent Directions algorithm.
%
% FD is based on Liberty et al.: https://arxiv.org/abs/1501.01711.pdf
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 17/04/2020
% 
% License: GPLv3
%
  fprintf('\n ** Running regular FD...\n');
  
  % initialisations
  m = 2 * ell;
  [~, cols] = size(Y);
  Br = zeros(m, cols);
  nz_row = 1;
  % set the projected data to be 0
  opt_out.Yr = 0;
  
  % check if we have parameters
  if nargin < 3
    params.no_final_err = 1;
    params.no_err = 1;
    params.use_blk_err = 1;
    params.err_blk_size = 100;
  else
    % check if no error flag is present
    if ~isfield(params, 'no_err')
      params.no_err = 1;
    end
    
    % check if use block error flag is present
    if ~isfield(params, 'use_blk_err')
      params.use_blk_err = 0;
    end
    
    % check if error block size is present
    if ~isfield(params, 'err_blk_size')
      params.err_blk_size = 100;
    end
    
    % check if no projected data flag is present
    if ~isfield(params, 'no_final_err')
      params.no_final_err = 0;
    end
  end

  % the number of rows
  numr = size(Y, 1);

  % initialise error metrics
  if params.no_err == 0
    % counter for the errors
    cnt = 1;
    % check if we use block errors or not
    if params.use_blk_err == 1
      opt_out.ErrFro = nan(1, floor(numr/params.err_blk_size));
      opt_out.T = nan(1, floor(numr/params.err_blk_size));
    else
      opt_out.ErrFro = nan(1, numr);
      opt_out.T = 1:numr;
    end
  % error tracking is disabled
  else
    opt_out.ErrFro = 0;
    opt_out.T = 0;
  end
  % start timing
  ts = tic;
  
  % loop through matrix
  for k = 1:numr
    % check if we need to squeeze
    if (nz_row >= m)
      % squeeze
      [Br, nz_row, ~] = fd_rotate_sketch(Br, ell);
    end
    % append the current values
    Br(nz_row, :) = Y(k, :);
    % increment the next zero row counter
    nz_row = nz_row + 1;

    % calcualte the error, if needed
    if params.no_err == 0
      if params.use_blk_err == 1
        if mod(k, params.err_blk_size) == 0
          y_c = Y(1:k, :);
          YrHat_c = y_c*(Br(1:ell, :)'*Br(1:ell, :));
          temp = sum(sum((y_c-YrHat_c).^2, 1));
          opt_out.ErrFro(cnt) = temp/k;
          opt_out.T(cnt) = k; cnt = cnt + 1;
        end
      else
        % calculate the reconstruction error
        y_c = Y(1:k, :);
        YrHat_c = y_c*(Br(1:ell, :)'*Br(1:ell, :));
        temp = sum(sum((y_c-YrHat_c).^2, 1));
        opt_out.ErrFro(k) = temp/k;
      end
    end
  end
  
  % set the final estimate of Yr, if we require the projected data.
  if params.no_err == 0
    if params.no_final_err == 0
      opt_out.Yr = Y*(Br(1:ell, :)'*Br(1:ell, :));
    end
  end
  
  % only return the subset of the sketch that is of value
  Bout = Br(1:ell, :);
  % calcualte the current trial execution delta
  opt_out.t = my_toc(ts);
end