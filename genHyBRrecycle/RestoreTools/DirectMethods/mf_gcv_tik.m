function alpha = mf_gcv_tik(s, bhat)
%GCV_TIK Choose GCV parameter for Tikhonov multiframe image deblurring.
%
%function alpha = mf_gcv_tik(s, bhat)
%
%         alpha = mf_gcv_tik(s, bhat);
%
%  This function uses generalized cross validation (GCV) to choose
%  a regularization parameter for Tikhonov filtering, when we are
%  using "Franklin's method.
%
%  Input:
%        s  Vector containing singular or spectral values
%             of the blurring matrix.
%     bhat  Vector containing the spectral coefficients of the blurred
%             image.
%
%  Output:
%    alpha  Regularization parameter.

% Reference: See Chapter 6, 
%            "Deblurring Images - Matrices, Spectra, and Filtering"
%            by P. C. Hansen, J. G. Nagy, and D. P. O'Leary,
%            SIAM, Philadelphia, 2006.

alpha = fminbnd(@mf_GCV, min(abs(s)), max(abs(s)), [], s, bhat);

  function G = mf_GCV(alpha, s, bhat)
    %
    %  This is a nested function that evaluates the GCV function for
    %  Tikhonov filtering.  It is called by fminbnd.
    %
    phi_d = 1 ./ (s + alpha^2);
    G = sum(abs(bhat.*phi_d).^2) / (sum(phi_d)^2);
  end
  
end