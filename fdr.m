function [Bout, alpha, ErrFro, T, Yr, t] = fdr(Y, ell, a_seed, no_err)
%FDR Find the robust frequent directions of a given matrix A in a
%streaming fashion
%
% Based on work from Liberty et al. and Luo et al.
%
% FD: https://arxiv.org/abs/1501.01711.pdf
% RFD: https://arxiv.org/pdf/1705.05067
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 30/12/2018
% 
% License: GPLv3
%
  fprintf('\n ** Running Robust FD...\n');

  % scope in global variables
  global use_blk_err

  % initialisations
  m = 2 * ell;
  [~, cols] = size(Y);
  Br = zeros(m, cols);
  nz_row = 1;
  
  % the number of rows and columns of Y
  [numr, numc] = size(Y);
  % initialise the identity matrix for the regulariser
  Id = eye(numc);

  % default block size
  blk_size = 100;
  cnt = 1;

  % default starting alpha value
  if nargin < 3
    a_seed = 0;
  end

  % set initial alpha
  alpha = 0;
  % previous alpha for the gradient
  alpha_prev = a_seed;

  % no error by default
  if nargin < 4
    no_err = 1;
  end

  % initialise error metrics
  if use_blk_err == 1
    ErrFro = nan(1, floor(numr/blk_size));
    T = nan(1, floor(numr/blk_size));
  else
    ErrFro = nan(1, numr);
    T = 1:numr;
  end

  % start timing
  ts = tic;
  
  % loop through matrix
  for k = 1:numr
    % check if we need to squeeze
    if (nz_row >= m)
      % squeeze
      [Br, nz_row, alpha] = fd_rotate_sketch(Br, ell, alpha_prev);
      % update the previous alpha regulariser value
      alpha_prev = alpha;
    end
    % append the current values
    Br(nz_row, :) = Y(k, :);
    % increment the next zero row counter
    nz_row = nz_row + 1;

    % calcualte the error, if needed
    if no_err == 0
      if use_blk_err == 1
        if mod(k, blk_size) == 0
          y_c = Y(1:k, :);
          YrHat_c = y_c*(Br(1:ell, :)'*Br(1:ell, :) + alpha*Id);
          temp = sum(sum((y_c-YrHat_c).^2, 1));
          ErrFro(cnt) = temp/k;
          T(cnt) = k; cnt = cnt + 1;
        end
      else
        % calculate the reconstruction error
        y_c = Y(1:k, :);
        YrHat_c = y_c*(Br(1:ell, :)'*Br(1:ell, :) + alpha*Id);
        temp = sum(sum((y_c-YrHat_c).^2, 1));
        ErrFro(k) = temp/k;
      end
    end

  end
  % also set the final estimate of Yr
  Yr = Y*(Br(1:ell, :)'*Br(1:ell, :) + alpha*Id);
  % only return the subset of the sketch that is of value
  Bout = Br(1:ell, :);
  % calcualte the current trial execution delta
  t = my_toc(ts);
end