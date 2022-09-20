function x = mldivide(A, b)
%
%  Overload backslash operation for mpsfMatrix0 object.
%
%  Solves Tikhonov regularization least squares problem,
%  where the regularization parameter is chosen via GCV.
%

%  J. Nagy  2/23/12

if ( isa(A, 'mpsfMatrix0') )
  PSFs = A.psf;
  %center = round((size(PSFs(:,:,1))+1)/2);
  center = A.center;
  [x, alpha] = mf_tik_fft(b, PSFs, center); 
else

  error('incorrect argument type')

end
