function [Bout, ErrFro, T, Yr, t] = fd(Y, ell, no_err)
%FD Find the frequent directions of a given matrix Y in a
%streaming fashion
%
% FD is based on Liberty et al.: https://arxiv.org/abs/1501.01711.pdf
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 30/12/2018
% 
% License: GPLv3
%
  fprintf('\n ** Running regular FD...\n');

  % scope in global variables
  global use_blk_err
  
  % initialisations
  m = 2 * ell;
  [~, cols] = size(Y);
  Br = zeros(m, cols);
  nz_row = 1;

  % default block size
  blk_size = 100;
  cnt = 1;

  % no error by default
  if nargin < 3
    no_err = 1;
  end

  % the number of rows
  numr = size(Y, 1);

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
      [Br, nz_row, ~] = fd_rotate_sketch(Br, ell);
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
          YrHat_c = y_c*(Br(1:ell, :)'*Br(1:ell, :));
          temp = sum(sum((y_c-YrHat_c).^2, 1));
          ErrFro(cnt) = temp/k;
          T(cnt) = k; cnt = cnt + 1;
        end
      else
        % calculate the reconstruction error
        y_c = Y(1:k, :);
        YrHat_c = y_c*(Br(1:ell, :)'*Br(1:ell, :));
        temp = sum(sum((y_c-YrHat_c).^2, 1));
        ErrFro(k) = temp/k;
      end
    end

  end
  % also set the final estimate of Yr
  Yr = Y*(Br(1:ell, :)'*Br(1:ell, :));
  % only return the subset of the sketch that is of value
  Bout = Br(1:ell, :);
  % calcualte the current trial execution delta
  t = my_toc(ts);
end