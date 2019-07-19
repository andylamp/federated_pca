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
function [W,  k, Proj, recon, errs, pcs, T, ErrFro] = SPIRIT(A, lambda, energy, k0, ...
  holdOffTime, silent, no_err)

if nargin < 7, no_err = 1; end
if nargin < 6, silent = 0; end
if nargin < 5, holdOffTime = 1; end
if nargin < 4, k0 = 3; end

% scope in global variables
global use_blk_err

% default block size
blk_size = 100;

n = size(A,2);
totalTime = size(A, 1);
Proj = zeros(totalTime,n); 
recon = zeros(totalTime,n);
pcs = zeros(n, 1);
%initialize w_i to unit vectors
W = eye(n);
d = 0.01*ones(n, 1);
m = k0; % number of eigencomponents

dec = 0;
inc = 0;
mmax = 0;

relErrors = zeros(totalTime, 1);

sumXSq=0;
sumYSq=0;

% initialise error metrics
if use_blk_err == 1
  ErrFro = nan(1, floor(totalTime/blk_size));
  T = nan(1, floor(totalTime/blk_size));
else
  ErrFro = nan(1, totalTime);
  T = 1:totalTime;
end

%incremental update W
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
  Y = W(:,1:m)' * A(t,:)'; %project to m-dimensional space
  xActual = A(t,:)'; %actual vector of the current time
  xProj = W(:,1:m) * Y; %reconstruction of the current time
  Proj(t,1:m) = Y; 
  recon(t,:) = xProj;
  xOrth = xActual - xProj;

  relErrors(t) = sum(xOrth.^2)/sum(xActual.^2);

  % update energy
  sumYSq = lambda * sumYSq + sum(Y.^2);
  sumXSq = lambda * sumXSq + sum(A(t,:).^2);
  
  % check the lower bound of energy level
  if(sumYSq < energy(1)*sumXSq && lastChangeAt < t - holdOffTime && m < n)
    lastChangeAt = t;
    m = m+1;
    % reinitialise the component
    c = zeros(1, n);
    c(m) = 1;
    W(:, m) = c;
    W(:,1:m) = grams(W(:,1:m));
    inc = inc + 1;
    if silent == 0
      fprintf('Increasing m to %d at time %d (ratio %6.2f)\n', m, t, 100*sumYSq/sumXSq);
    end
  % check the upper bound of energy level
  elseif (sumYSq > energy(2)*sumXSq && lastChangeAt < t - holdOffTime && m < n && m>1)
    lastChangeAt = t;
    m = m-1;
    dec = dec + 1;
    if silent == 0
      fprintf('Decreasing m to %d at time %d (ratio %6.2f)\n', m, t, 100*sumYSq/sumXSq);  
    end
  end
  
  if(mmax < m)
      mmax = m;
  end
  % update the pcs
  pcs(t) = m;
  
    % this is disabled in the speed run
  if no_err == 0
    % check if we use block error
    if use_blk_err == 1        
      if mod(t, blk_size) == 0
        y_c = A(1:t, :)';
        YrHat_c = (W*W')*y_c;
        temp = sum(sum((y_c-YrHat_c).^2, 1));
        ErrFro(cnt) = temp/t;
        T(cnt) = t; cnt = cnt + 1;
      end
    else
      % Frobenius norm incremental error, per column in the block at 
      % (1:k*B) normalised with current T.
      y_c = A(1:t, :)';
      YrHat_c = (W*W')*y_c;
      temp = sum(sum((y_c-YrHat_c).^2, 1));
      ErrFro(t) = temp/t;
    end
  end
    
end

% check if we are in a speed run where errors
% are disabled.
if no_err == 0
  ErrFro = ErrFro(:);
  T = T(:);
  if use_blk_err == 1
    % avoid the first block error
    T = [0; T];
    ErrFro = [0; ErrFro];
  end
else
  ErrFro = 0;
  T = 0;
end

fprintf(' -- SPIRIT: Final m: %d at time %d (ratio %6.2f) [total decs: %d, total incs: %d, max m: %d]\n', m, t, 100*sumYSq/sumXSq, dec, inc, mmax);  

% set outputs
W(:,1:m) = grams(W(:,1:m));
W = W(:,1:m);
k = m;
errs = relErrors;

%figure;
%plot(dtrack);
%plot(Proj(1:30000,:));
%plot(Proj(1:7000,:));
%plot(Proj);
%title('SPIRIT Projections');
% figure;
% plot(errs);

