function [Um, Sm] = fpca_subspace_merge(U1, S1, U2, S2, lambda1, lambda2, r, type)
%FPCA_SUBSPACE_MERGE Merges two subspaces of equal or different ranks.
% Requires both (orthogonal) subspaces along with their singular values. 
% Futher, optionally, we can provide weighting factors for each subspace 
% (lambdas {1, 2}), the desired rank as well as the type of 
% merge algorithm used.
%
% This method currently supports three merging algorithm types -- by  
% default we use type 3.
%
% Supported types:
%
% 1) purely SVD based (good if d and r is small)
%
% 2) slightly more efficient QR + SVD based (good for medium d and r)
% 
% 3) best memory effiency with excellent preservation of Singular 
%    Values using again a combination or QR + SVD but the block matrix 
%    that the SVD is applied upon is munch smaller.
%
% Based on work of Grammenos et al.: https://arxiv.org/abs/1907.08059
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 26/04/2020
% 
% License: GPLv3
%
  % find the limits
  r1 = size(U1, 2);
  r2 = size(U2, 2);
  
  % check if S1 is a vector as we need to expand
  if isvector(S1)
    S1 = diag(S1);
  end
  
  % check if S2 is a vector as we need to expand
  if isvector(S2)
    S2 = diag(S2);
  end
  
  % ensure we have r x r singular values
  S1 = S1(1:r1, 1:r1);
  S2 = S2(1:r2, 1:r2);
  rmax = max(r1, r2);
  
  % target type of algo - default we use type 3
  if nargin < 8
    type = 3;
  end
  
  % target rank, otherwise default to max of the two
  if nargin < 7
    r = rmax;
  elseif r > rmax 
    r = rmax;
  end

  % check if our forgetting factors are present, otherwise
  % set to default of 1
  if nargin < 5
    lambda1 = 1;
    lambda2 = 1;
  end
  
  % check if we have nan values
  if isnan(U1)
    Um = U2;
    Sm = S2;
    return;
  elseif isnan(U2)
    Um = U1;
    Sm = S1;
    return
  end
  
  % ts = tic;
  % "optimal"  
  if type == 1
    [Um, Sm, ~] = svds([U1*S1, U2*S2], 2*r);
  % somewhat efficient
  elseif type == 2
    [Qp, Rp] = qr([lambda1*(U1*S1), lambda2*U2*S2], 0);
    [Ur, Sr, ~] = svds(Rp, r);
    qq = Qp(:, 1:min(r1+r2, rmax));
    uu = Ur(1:(min(r1+r2, rmax)), 1:r);
    Um = qq*uu;
    Sm = Sr;
  % naive
  elseif type == 3
    % get the multiplication of the subspace and descent to 
    % r1 x r2 space.
    z_k = U1'*U2;
    % perform economy QR to get the U' and R
    [Up, Rp] = qr(U2-U1*z_k, 0);
    % construct the block matrix of (r1+r2) x (r1+r2) space
    blk_mat_s = [lambda1*S1, z_k*S2; zeros(r2, r1), lambda2*Rp*S2];
    % get the svds from the block mat using the target rank
    % which by default is the max(r1, r2)
    [Ur, Sr, ~] = svds(blk_mat_s, 2*r);
    % get the reference R1 & R2 from Ur
    R1 = Ur(1:r1, 1:r);
    R2 = Ur(1:r2, 1:r);
    % finally set the outputs
    Um = U1*R1 + Up*R2;
    Sm = Sr(1:r, 1:r);
  % unrecognised, signal the error
  else
    error('unknown type, types have to be between 1-3');
  end
  %my_toc(ts);
  Um = Um(:, 1:r);
  Sm = Sm(1:r, 1:r);
end

