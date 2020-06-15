function [beta] = beta_xn(x, n, e_p, delta)
%BETA_XN Used to compute the xn values for the test d/B ratio.
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 03/06/2020
% 
% License: GPLv3
%
  % first part
  p1 = ((x * n)+1) / (n * e_p);
  % the logarithm inside.
  log_in = ((x * n) * (x * n + 1)) / (2 * delta * sqrt(2 * pi));
  % the second part.
  p2 = 2 * log(log_in);
  % the third part.
  p3 = 1 / (n * sqrt(e_p));
  % add them all together.
  beta = (p1 * p2) + p3;
end

