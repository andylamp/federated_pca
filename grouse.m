function [U, R, err_reg, T, ErrFro] = grouse(I, J, S, numr, numc, params) 
%GROUSE Grassman Rank-One Update Subspace Estimation matrix completion code
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  GROUSE (Grassman Rank-One Update Subspace Estimation) matrix completion code 
%  by Ben Recht and Laura Balzano, February 2010.
%
%  Given a sampling of entries of a matrix X, try to construct matrices U
%  and R such that U is unitary and UR' approximates X.  This code 
%  implements a stochastic gradient descent on the set of subspaces.
%
%  Inputs:
%       (I,J,S) index the known entries across the entire data set X. So we
%       know that for all k, the true value of X(I(k),J(k)) = S(k)
%
%       numr = number of rows
%       numc = number of columns
%           NOTE: you should make sure that numr<numc.  Otherwise, use the
%           transpose of X
%       
%       max_rank = your guess for the rank
%
%       step_size = the constant for stochastic gradient descent step size
%
%       maxCycles = number of passes over the data
%
%       Uinit = an initial guess for the column space U (optional)
%
%   Outputs:
%       U and R such that UR' approximates X.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% assign wrapper parameters
err_blk_size = params.err_blk_size;
max_rank = params.max_rank;
step_size = params.step_size;
max_cycles = params.max_cycles;
en_sparse = params.en_sparse;
rperm = params.rperm;
no_err = params.no_err;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matlab specific data pre-processing
%

% Form some sparse matrices for easier matlab indexing
%values = sparse(I,J,S,numr,numc);
if en_sparse == 1
    values = sparse(I,J,S,numr,numc);
    Indicator = sparse(I,J,1,numr,numc);
else
    values = S;
    Indicator = ones(numr, numc);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Main Algorithm
%

if ~isfield(params, 'Uinit')
  % initialize U to a random r-dimensional subspace 
  U = orth(randn(numr,max_rank)); 
else
  U = params.Uinit;
end

err_reg = zeros(max_cycles*numc,1);

fprintf("\n\t ** Max Passes: %d", max_cycles);
if rperm == 1
  fprintf('\n\t ** Random Column Permutation is ENABLED\n\n');
else
  fprintf('\n\t ** Random Column Permutation is DISABLED\n\n');
end

% check if need to pre-allocate the error metrics
if no_err == 0
  % initialise error metrics
  if params.use_blk_err == 1
    ErrFro = nan(max_cycles, floor(numc/err_blk_size));
    T = nan(max_cycles, floor(numc/err_blk_size));
    % counter for block error (if enabled)
    cnt = 1;
  else
    ErrFro = nan(max_cycles, numc);
    T = 1:max_cycles*numc;
  end
end

for outiter = 1:max_cycles
    
    fprintf('\t !! Pass %d...\n', outiter);
    
    % create a random ordering of the columns for the current pass over the
    % data.
    if rperm == 1
        col_order = randperm(numc);
        Yc = values(:, col_order);
    else 
        % no random permutation
        col_order = 1:numc;
        Yc = values;
    end
    
for k=1:numc
      
    % Pull out the relevant indices and revealed entries for this column
    idx = find(Indicator(:, col_order(k)));
    v_Omega = values(idx, col_order(k));
    U_Omega = U(idx,:);    

    
    % Predict the best approximation of v_Omega by u_Omega.  
    % That is, find weights to minimize ||U_Omega*weights-v_Omega||^2
    
    weights = U_Omega\v_Omega;
    norm_weights = norm(weights);
    
    % Compute the residual not predicted by the current estimate of U.

    residual = v_Omega - U_Omega*weights;       
    norm_residual = norm(residual);
    
    % This step-size rule is given by combining Edelman's geodesic
    % projection algorithm with a diminishing step-size rule from SGD.  A
    % different step size rule could suffice here...        
    
    sG = norm_residual*norm_weights;
    err_reg((outiter-1)*numc + k) = norm_residual/norm(v_Omega);
    t = step_size*sG/( (outiter-1)*numc + k );
   
    % Take the gradient step.    
    if t<pi/2 % drop big steps        
        alpha = (cos(t)-1)/norm_weights^2;
        beta = sin(t)/sG;

        step = U*(alpha*weights);
        step(idx) = step(idx) + beta*residual;

        U = U + step*weights';
    end
    % this is disabled in the speed run
    if no_err == 0
      % check if we use block error
      if params.use_blk_err == 1        
        if mod(k, err_blk_size) == 0
          y_c = Yc(:, 1:k);
          YrHat_c = (U*U')*y_c;
          temp = sum(sum((y_c-YrHat_c).^2, 1));
          ErrFro(outiter, cnt) = temp/k;
          T(outiter, cnt) = k; cnt = cnt + 1;
        end
      else
        % Frobenius norm incremental error, per column in the block at 
        % (1:k*B) normalised with current T.
        y_c = Yc(:, 1:k);
        YrHat_c = (U*U')*y_c;
        temp = sum(sum((y_c-YrHat_c).^2, 1));
        ErrFro(outiter, k) = temp/k;
      end
    end
end

end


% check if we are in a speed run where errors
% are disabled.
if no_err == 0
  ErrFro = ErrFro(:);
  T = T(:);
  if params.use_blk_err == 1
    % avoid the first block error
    T = [0; T];
    ErrFro = [0; ErrFro];
  end
else
  ErrFro = 0;
  T = 0;
end

% Once we have settled on our column space, a single pass over the data
% suffices to compute the weights associated with each column.  You only
% need to compute these weights if you want to make predictions about these
% columns.
%fprintf('Find column weights...');
R = zeros(numc, max_rank);
for k= 1:numc     
    % Pull out the relevant indices and revealed entries for this column
    idx = find(Indicator(:, k));
    v_Omega = values(idx, k);
    U_Omega = U(idx, :);
    % solve a simple least squares problem to populate R
    R(k, :) = (U_Omega\v_Omega)';
end

