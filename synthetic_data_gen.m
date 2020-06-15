function [Y, spectrum] = synthetic_data_gen(n, T, params)
%SYNTHETIC_DATA_GEN Generates T random vectors in R^{1 x n}.
% The vectors are drawn from Power Law distribution with parameters 
% `lambda` and `alpha` which indicate the spectrum of the singular values.
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 12/06/2020
% 
% License: GPLv3
% 

  % now decide which distributions we have and generate their respective 
  % spectrum, starting with Power Law (pl).
  if params.spectrum_type == "pl"
    % check if lambda is a parameter
    if ~isfield(params, "lambda")
      params.lambda = 1;
    end
    
    % check if alpha is a parameter
    if ~isfield(params, "alpha")
      params.alpha = 1;
    end
    stitle = sprintf("Power law scree plot of alpha: %6.2f", params.alpha);
    spectrum = params.lambda*(1:n).^(-params.alpha);
  % The other distribution we have is singular - i.e.: s1=1, then others 
  % are a sufficiently small value such as 0.01.
  elseif params.spectrum_type == "singular"
    % check if lambda is a parameter
    if ~isfield(params, "lambda")
      params.lambda = 0.01;
    end

    stitle = sprintf("Singular scree plot (\lambda: %3.2f)", params.lambda);
    spectrum = [1, params.lambda*ones(1, n-1)];
  % finally this is the spectrum used in mod-sulq paper and is used only
  % for testing.
  elseif params.spectrum_type == "mod-sulq"
    stitle = sprintf("MOD-SuLQ paper scree plot");
    spectrum = [0.5, 0.3, 0.04, 0.03, 0.02, 0.01, 0.004, 0.003, 0.001, 0.001];
    % check if n is valid
    if n ~= size(spectrum, 2)
      error("spectrum dimensions must match with ambient");
    end
  end
  
  % check for silent flag, default 1
  if ~isfield(params, "silent")
    params.silent = 1;
  end
  
  % check if we have a rand type field
  if ~isfield(params, "rand_type")
    params.rand_type = "normal";
  end

  % generate Sigma
  Sigma = diag(spectrum);
  if params.silent == 0
    figure; plot(1:n, diag(Sigma).^2); title(stitle);
  end
  
  % generate the synthetic data based on the selected spectrum type
  if params.rand_type == "uniform"
    % get the basis
    [U, ~] = qr(rand(n));
    % given S and Sigma generate the dataset (Y)
    Y = (U * Sigma * rand(n, T)) / sqrt(T-1);
  elseif params.rand_type == "normal"
    % get the basis
    [U, ~] = qr(randn(n));
    % U = orth(randn(n));
    % given S and Sigma generate the dataset (Y)
    Y = (U * Sigma * randn(n, T)) / sqrt(T-1);
  end
  % find the max norm (to have a unit while preserving the distribution)
  yn = max(sqrt(sum(Y.^2)));
  % now divide the Y with that value
  Y = Y./yn;
end

