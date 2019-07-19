function [T, ErrFro, U, YrHat, t] = mitliag_pm(Y, r, blkSize, floor_mul, no_err)
%%MITLIAG_PM This is an implementation of the Power Method as described in  
% "Memory Limited, Streaming PCA" (https://arxiv.org/pdf/1307.0032.pdf)
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 30/12/2018
% 
% License: GPLv3
%
  fprintf('\n ** Running Mitliagkas Power Method...\n')
  
  % scope in global variables
  global use_blk_err
  
  % ambient dim defined by our input matrix Y
  [n, Ti] = size(Y);
  
  % check we calculate the error (disabled for speed runs)
  if nargin < 5
    no_err = 0;
  end

  % check if we have a floor multiplier as an argument
  if nargin < 4
    floor_mul = 2;
  end

  % check if n < r or n == 1
  if n == 1 || n < r
    error(" ** ERR: Ambient dimension must be > 1 and r < n **");
  end

  % Random initialization
  SrHat = orth(randn(n));
  SrHatTemp = SrHat(:, 1:r);

  % power method configuration
  % set the block, depending on argument
  if nargin < 3
    %B = 2*n;
    B = n;
  else
    if blkSize < n
      fprintf("\n !! WARN: Block size must be at least n !!\n");
      B = n;
    else
      B = blkSize;
    end
  end
  
  % check if Ti < B, in which case we cannot run it
  if Ti < B
    error("\n Block size must be lower than the number of columns");
  end
  
  K = floor(Ti/B);              % Number of blocks
  cnt = 1;                      % index for error block align 
  
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
  
  % output the number of blocks and their size
  fprintf([' ** Total number of blocks (k): %d ', ...
    'with Block size of: %d\n'], K, B);
  
  % start timing
  ts = tic;
  for k = 1:K
    % update the previous subspace estimation
    SrHat_old = SrHatTemp;
    % indices for current block
    min_t = ((k-1)*B)+1;
    max_t = k*B;
    % initialise the temporary place-holder
    temp = zeros(n,r);
    for i = 1:B
      index = (k-1)*B+i;
      temp = temp+1/B*(Y(:,index)*(Y(:,index))')*SrHat(:, 1:r);
    end 
    [SrHat, ~, ~] = svds(temp, floor(floor_mul * r));
    SrHatTemp = SrHat(:, 1:r);

    % Mitliagkas Power Method Errors start

    if no_err == 0
      % use the block error
      if use_blk_err == 1      
        % generate current Yr estimation
        YrHat_c = (SrHat*SrHat')*Y(:, 1:k*B);
        % Frobenius norm incremental error, per block at 
        % at (k-1)B:kB normalised with current T.
        temp = sum(sum((Y(:, 1:max_t)-YrHat_c).^2, 1));
        ErrFro(cnt) = temp/max_t;
        T(cnt) = max_t;
        cnt = cnt + 1;
      else
        % generate current Yr estimation
        YrHat_c = (SrHat_old*SrHat_old')*Y(:, 1:k*B);
        % Frobenius norm incremental error, per column in the block  
        % at (k-1)B:kB normalised with current T.
        temp = cumsum(sum((Y(:, 1:max_t)-YrHat_c).^2, 1));
        ErrFro(cnt:cnt+B-1) = temp(min_t:max_t)./(min_t:max_t);
        cnt = cnt + B;
      end
    end

    % Mitliagkas Power Method errors end
    
  end
  % check what to do in case of error generation
  if no_err == 0
    YrHat = YrHat_c;  % current YrHat estimation
  else
    YrHat = (SrHatTemp*SrHatTemp')*Y(:, 1:k*B);
  end
  % update the other estimates
  U = SrHatTemp;   % current U estimation
  t = my_toc(ts);  % current elapsed time 
end