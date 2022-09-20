function [X, alpha] = mf_tik_fft(B, PSF, center, alpha)
%TIK_FFT Tikhonov multiframe image deblurring using the FFT algorithm.
%
%function [X, alpha] = mf_tik_fft(B, PSF, center, alpha)
%
%            X = mf_tik_fft(B, PSF, center);
%            X = mf_tik_fft(B, PSF, center, alpha);
%   [X, alpha] = mf_tik_fft(B, PSF, ...);
%
%  Compute multiframe restoration using an FFT-based Tikhonov filter, 
%  with the identity matrix as the regularization operator.
%
%  Input:
%        B  3D Array containing blurred images.
%      PSF  3D Array containing PSFs; same size as B.
%   center  [row, col] = indices of center of PSFs.
%    alpha  Regularization parameter.
%             Default parameter chosen by generalized cross validation.
%
%  Output:
%        X  Array containing computed restoration.
%    alpha  Regularization parameter used to construct restoration.

% Reference: See Chapter 6, 
%            "Deblurring Images - Matrices, Spectra, and Filtering"
%            by P. C. Hansen, J. G. Nagy, and D. P. O'Leary,
%            SIAM, Philadelphia, 2006.

%
% Check number of inputs and set default parameters.
%
if (nargin < 3)
   error('B, PSF, and center must be given.')
end
if (nargin < 4)
   alpha = [];
end

%
% Use the FFT to compute the eigenvalues of the BCCB blurring matrices,
% and sum.  Also, accumulate right hand side.
%
% In the multi-frame case, we assume all PSFs have the same center.
%
S = zeros(size(PSF,1), size(PSF,2));
bhat = zeros(size(B,1), size(B,2));
Sfft = fft2( circshift(PSF, 1-center) );
S = sum(abs(Sfft).^2,3);
bhat = sum(conj(Sfft).*fft2(B),3);
s = S(:);

%
% If a regularization parameter is not given, use GCV to find one.
% Note at this point we have alread formed the normal equations,
% so we should use "Franklin's" method. (See P.C. Hansen book).
%
bhat = bhat(:);
if (ischar(alpha) || isempty(alpha))
  alpha = mf_gcv_tik(s, bhat);
  %
  %  It seems GCV is choosing a parameter too small, but if we use
  %  sqrt of the chosen value, it works a bit better.
  %
  alpha = sqrt(alpha);
end 
  
%
% Compute the Tikhonov regularized solution, using Franklin's method.
%
D = s + abs(alpha)^2;
xhat = bhat ./ D;
xhat = reshape(xhat, [size(B,1), size(B,2)]);
X = real(ifft2(xhat));