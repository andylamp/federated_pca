function [Ur, Sr, V] = ssvdr_block_update(B, r, Up, Sp, V)
%SSVDR_BLOCK_UPDATE Update the estimate of r-truncated SVD. 
%
% Based on work of Grammenos et al.: https://arxiv.org/abs/1907.08059
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 19/04/2020
% 
% License: GPLv3
%
  % check if its the initial block
  if nargin < 4  
    Up = NaN;
    Sp = NaN;
  end
  
  % check if we request vt
  if nargin < 5
    V = NaN;
  end
  
  % check if we return Vt
  if ~isnan(V) 
    ret_vt = 1;
  else
    ret_vt = 0;
  end
  
  % check if we have the first block
  if isnan(Up)
    if ret_vt == 1
      [Ur, Sr, V] = svds(B, r);
    else
      [Ur, Sr, ~] = svds(B, r);
    end
  % we do not, process the next block
  else
    % construct the q_k (projections)
    q_k = Up(:, 1:r)'*B;
    % reconstruction energy
    cur_recon = Up(:, 1:r)*q_k;
    % construct the z_k
    z_k = B - cur_recon;

    % get the (economy) QR of z_k
    [s_k, v_k] = qr(z_k, 0);

    % now construct the following block matrix as is shown in the
    % algorithm in our paper:
    %
    %           |      G_k      q_k |
    % blk_mat = |                   |
    %           | zeros(zr, r)  v_k |
    %
    zr = size(v_k, 1);
    %zr = min(b, size(v_k, 1));
    % grab the previous diagonal eigenvalues, while check if it is a vector
    if isvector(Sp)
      % expand into a proper singular value matrix
      G_k = diag(Sp);
    else
      % it is already expanded.
      G_k = Sp;
    end
    %
    %cc_top = [G_k(1:r, 1:r), q_k(:, :)]; %size(cc_top)
    %cc_bot = [zeros(zr, r), v_k]; %size(cc_bot)
    %
    %k
    blk_mat_k = [ G_k(1:r, 1:r), q_k(:, :); zeros(zr, r), v_k ];

    % now take the r-svds of that matrix
    if ret_vt == 1
      [u_k, Sr, vt] = svds(blk_mat_k, 2*r);
    else
      [u_k, Sr, ~] = svds(blk_mat_k, 2*r);
    end
    
    % now update the actual Ur estimation
    blk_st = [Up(:, 1:r), s_k];
    u_k_sz = size(u_k, 1);
    Ur = blk_st*u_k(1:u_k_sz, 1:r);
    
    % check if we require the projected data
    if ret_vt == 1
      bs = size(B, 2);
      V = [V, zeros(size(V, 1), bs); zeros(bs, r), eye(bs)]*vt;
    end
    
  end
end

