function [O] = gm_orth(A)
%GM_ORTH Vectorised Gram-Schmidt orthonormoalisation.
%
% In order for this to work the columns of A are assumed to be indeed
% linearly indepedent; otherwise, perhaps try to use a more stable
% algorithms such as qr.
%
%
% Author Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date 06/05/2020
% 
% License: GPLv3
% 
  % sanity checks.
	if nargin~=1 || size(A,2)<2 
    error(" !! ERR: Input matrix needs to have at least 2 columns."); 
  end
  % preallocate the variables for speed
	O = 0*A;
  Z = 0*A;
  % perform the vectorised orthonormalisation
	for i=1:size(A,2)
		Z(:,i) = A(:,i) - sum(O(:,1:i) * diag(A(:,i)' * O(:,1:i)), 2);
		O(:,i)= Z(:,i) / norm(Z(:,i));
	end	
end
