function [B_out, nz_row, alpha] = fd_rotate_sketch(B_in, ell, alpha_prev)
%FD_ROTATE_SKETCH the main shirk and rotation that's performed when
%receiving a new row and our buffer is full from our stream.
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 30/12/2018
% 
% License: GPLv3
%
    % initialise output
    B_out = B_in;
    % default alpha value is zero, if we don't use it
    if nargin < 3
      alpha_prev = 0;
    end
    % calculate the svds of B
    [~, S, Vt] = svd(B_in);
    Sd = diag(S);
    [sd_rows, ~] = size(Sd);
    if sd_rows >= ell
        % take the square error or the last row compared to the sketch
        shrunk_sketch = sqrt(Sd(1:ell, :).^2 - Sd(ell).^2);
        % update the sketch
        B_out(1:ell, :) = diag(shrunk_sketch) * Vt(1:ell, :);
        % zero out the last row
        B_out(ell + 1, :) = 0;
        nz_row = ell + 1;
    else
        % update the portion of the sketch
        B_out(1:sd_rows, :) = S * Vt(1:sd_rows, :);
        % zero out the last row of the sketch
        B_out(1:sd_rows + 1, :) = 0;
        nz_row = sd_rows + 1;
    end
    % calculate the new regulariser vector
    alpha = alpha_prev + (Sd(ell)^2)/2;
end