%% Function to construct a Gaussian blur
function [PSF, center] = Gauss( gamma , N )
%
%     [PSF, center] = Gauss( gamma , N )
%
% This function constructs the Gaussian blur PSF, and the center of the PSF
%
%   Input:
%        gamma - parameters defining the blur
%            N - size fo the PSF
%
%    Output:
%          PSF - Gaussian point spread function
%      center  - center of the PSF
%
%  J. Chung 1/24/09
%

m = N;
n = N;

%
% Set up grid points to evaluate the Gaussian function.
%
x = -fix(n/2):ceil(n/2)-1;
y = -fix(m/2):ceil(m/2)-1;
[X,Y] = meshgrid(x,y);

%
% Compute the Gaussian PSF
%
if length(gamma) == 1
  s1 = gamma(1);
  PSF = exp( -(x.^2)/(2*s1^2));
  PSFsum = sum(PSF(:));
elseif length(gamma) == 2
  s1 = gamma(1); s2 = gamma(2);
  PSF = exp( -(X.^2)/(2*s1^2) - (Y.^2)/(2*s2^2) );
  PSFsum = sum(PSF(:));
  
elseif length(gamma)==3
  s1 = gamma(1); s2 = gamma(2); s3 = gamma(3);
  num = -((X.^2)*(s1^2) + (Y.^2)*(s2^2) - 2*(s3^2)*(X.*Y));
  den = 2*(s1^2 * s2^2 - s3^4);
  PSF = exp( num / den );
  PSFsum = sum(PSF(:));
  
end
  
PSF = PSF/PSFsum;
center = [m/2+1,n/2+1];
end