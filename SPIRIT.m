%Simulation on updating w (works better for data with mean-zero)
%function [W, k, Proj, recon] = SPIRIT(A, lambda, energy, k0, holdOffTime)
%dynamic maintain the number of output units according to the energy level
%Inputs:
%A: data matrix
%-lambda: forgetting factor between [0,1], 1 means putting equal weight on
%         the past, otherwise the past is exponentially decayed by lambda
%-energy: an interval within  [0,1], which specifies the upper and lower bounds of energy level
%-k0: initial number of hidden variables, default=3
%-holdOffTime: the minimal timestamps before next change on the
%              number of hidden variables
%
%Output: 
%-W: eigvector matrix such as 1st eigvector W(:,1), kth eigvector W(:,k)
%-k: the max number of eigenvectors,
%-Proj: the lowe-dim projections;
%-recon: the reconstructions of A; 
%Example:
% >> A = sin((1:1000)/20)'*rand(1,10);
% >> plot(A)
% >> A = [sin((1:1000)/10) sin((1:1000)/50)]'*rand(1,10);
% >> plot(A)
% >> [W,k,Proj,recon] = SPIRIT(A,1,[0.95,0.98]);
% Decreasing m to 2 at time 12 (ratio 205.84)
% Decreasing m to 1 at time 23 (ratio 175.68)
% >> plot(recon)
%----------------------------------------------------------------
% Copyright: 2006,
%         Spiros Papadimitriou, Jimeng Sun, Christos Faloutsos.
%         All rights reserved.
% Please address questions or comments to: jimeng@cs.cmu.edu
%----------------------------------------------------------------
function [W, opt_out] = SPIRIT(A, lambda, energy, params)
%SPIRIT An implementation of SPIRIT, Papadimitriou et al.
%

  % if we do not have params, put everything in default.
  if nargin < 4
    params.no_err = 1;
    params.verbose = 0;
    params.holdoff_time = 1;
    params.use_blk_err = 0;
    params.err_blk_size = 100;
    params.k0 = 3;
    params.no_final_err = 1;
    params.adaptive = 1;
    params.use_qr = 0;
  % params was supplied, check if we have all entries otherwise
  % set missing to defaults.
  else
    % check if we have to use qr (necessary for robustness but lower perf.)
    if ~isfield(params, 'use_qr')
      params.use_qr = 0;
    end
    
    % check if adaptive is present
    if ~isfield(params, 'adaptive')
      params.adaptive = 1;
    end
    
    % check if no error flag is present
    if ~isfield(params, 'no_err')
      params.no_err = 1;
    end
    
    % check if verbose flag is present
    if ~isfield(params, 'verbose')
      params.verbose = 0;
    end
    
    % check if hold off time is present
    if ~isfield(params, 'holdoff_time')
      params.holdoff_time = 1;
    end
    
    % check if use block error flag is present
    if ~isfield(params, 'use_blk_err')
      params.use_blk_err = 0;
    end
    
    % check if an initial rank is provided
    if ~isfield(params, 'k0')
      params.k0 = 3;
    end
    
    % check if the error block size is set
    if ~isfield(params, 'err_blk_size')
      params.err_blk_size = 100;
    end
    
    % check if no projected data flag is present
    if ~isfield(params, 'no_final_err')
      params.no_final_err = 1;
    end
  end

n = size(A,2);
totalTime = size(A, 1);

ts = tic;

%initialize w_i to unit vectors
W = eye(n);
d = 0.01*ones(n, 1);
% number of eigencomponents
m = params.k0; 

dec = 0;
inc = 0;
mmax = 0;

sumXSq=0;
sumYSq=0;

opt_out.Yr = 0;

% initialise error metrics
if params.no_err == 0
  % basic errors
  if params.use_blk_err == 1
    opt_out.ErrFro = nan(1, floor(totalTime/params.err_blk_size));
    opt_out.T = nan(1, floor(totalTime/params.err_blk_size));
    cnt = 1;
  else
    opt_out.ErrFro = nan(1, totalTime);
    opt_out.T = 1:totalTime;
  end
  % misc errors
  opt_out.Proj = zeros(totalTime,n); 
  opt_out.recon = zeros(totalTime,n);
  opt_out.relErrors = zeros(totalTime, 1);
  opt_out.pcs = zeros(n, 1);
else
  opt_out.ErrFro = 0;
  opt_out.T = 0;
  opt_out.Proj = 0; 
  opt_out.recon = 0;
  opt_out.pcs = 0; 
end
  
% incremental update W
lastChangeAt = 1;

for t = 1:totalTime
  % update W for each y_t
  x = A(t,:)';
  for j = 1:m
     [W(:,j), d(j), x] = updateW(x, W(:,j), d(j), lambda);
     %Wj = W(:,j);
  end
    
  %W(:,1:m) = grams(W(:,1:m));
  %compute low-D projection, reconstruction and relative error
  Y = W(:, 1:m)' * A(t,:)'; %project to m-dimensional space
  
  % only compute the errors if we have to
  if params.no_err == 0
    xActual = A(t,:)'; %actual vector of the current time
    xProj = W(:,1:m) * Y; %reconstruction of the current time
    opt_out.Proj(t,1:m) = Y; 
    opt_out.recon(t,:) = xProj;
    xOrth = xActual - xProj;

    opt_out.relErrors(t) = sum(xOrth.^2)/sum(xActual.^2);
  end

  % update energy
  sumYSq = lambda * sumYSq + sum(Y.^2);
  sumXSq = lambda * sumXSq + sum(A(t,:).^2);
  
  % check the lower bound of energy level
  if params.adaptive == 1
    if(sumYSq < energy(1)*sumXSq && lastChangeAt < t - params.holdoff_time && m < n)
      lastChangeAt = t;
      m = m+1;
      % reinitialise the component
      c = zeros(1, n);
      c(m) = 1;
      W(:, m) = c;
      
      % orthnormalise
      if params.use_qr == 1
        W(:, 1:m) = qr(W(:, 1:m), 0);
      else
        W(:, 1:m) = grams(W(:, 1:m));
      end
      
      inc = inc + 1;
      if params.verbose == 1
        fprintf('Increasing m to %d at time %d (ratio %6.2f)\n', m, t, 100*sumYSq/sumXSq);
      end
    % check the upper bound of energy level
    elseif (sumYSq > energy(2)*sumXSq && lastChangeAt < t - params.holdoff_time && m < n && m>1)
      lastChangeAt = t;
      m = m-1;
      dec = dec + 1;
      if params.verbose == 1
        fprintf('Decreasing m to %d at time %d (ratio %6.2f)\n', m, t, 100*sumYSq/sumXSq);  
      end
    end
  end
  
  if(mmax < m)
      mmax = m;
  end
  % update the pcs
  opt_out.pcs(t) = m;
  
    % this is disabled in the speed run
  if params.no_err == 0
    % check if we use block error
    if params.use_blk_err == 1        
      if mod(t, params.err_blk_size) == 0
        y_c = A(1:t, :)';
        YrHat_c = (W(:, 1:m)*W(:, 1:m)')*y_c;
        temp = sum(sum((y_c-YrHat_c).^2, 1));
        opt_out.ErrFro(cnt) = temp/t;
        opt_out.T(cnt) = t; cnt = cnt + 1;
      end
    else
      % Frobenius norm incremental error, per column in the block at 
      % (1:k*B) normalised with current T.
      y_c = A(1:t, :)';
      YrHat_c = (W(:, 1:m)*W(:, 1:m)')*y_c;
      temp = sum(sum((y_c-YrHat_c).^2, 1));
      opt_out.ErrFro(t) = temp/t;
    end
  end
    
end

% check if we are in a speed run where errors
% are disabled.
if params.no_err == 0
  opt_out.ErrFro = opt_out.ErrFro(:);
  opt_out.T = opt_out.T(:);
  if params.use_blk_err == 1
    % avoid the first block error
    opt_out.T = [0; opt_out.T];
    opt_out.ErrFro = [0; opt_out.ErrFro];
  end
end

fprintf(' -- SPIRIT: Final m: %d at time %d (ratio %6.2f) [total decs: %d, total incs: %d, max m: %d]\n', m, t, 100*sumYSq/sumXSq, dec, inc, mmax);  

% check if we require the projected data
if params.no_final_err == 0
    opt_out.Yr = (W(:, 1:m)*W(:, 1:m)')*A;
end

% set outputs
if params.use_qr == 1
  W = qr(W(:, 1:m), 0);
else
  W = grams(W(:, 1:m));
end
% set the final rank
opt_out.final_rank = m;
% set the final time
opt_out.t = my_toc(ts);

