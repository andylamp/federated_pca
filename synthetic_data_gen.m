function [ Y, T, Sigma ] = synthetic_data_gen(n, T, lambda, alpha, silent)
%SYNTHETIC_DATA_GEN function that generates T random vectors in R^{1 x n} 
% drawn from Power Law distribution with parameters `lambda` and `alpha`
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 30/12/2018
% 
% License: GPLv3
% 

  % check for silent flag, default 1
  if nargin < 5
    silent = 1;
  end

  % generate the singular spectrum
  dd = lambda*(1:n).^(-alpha);
  % generate Sigma
  Sigma = diag(dd);
  if silent == 0
    stitle = sprintf("Scree plot of alpha: %6.2f", alpha);
    figure; plot(1:n, diag(Sigma).^2); title(stitle);
  end
  % random initialization of S basis
  S = orth(randn(n));
  % given S and Sigma generate the dataset (Y)
  Y = (S * Sigma * randn(n, T))/sqrt(T-1);

end

