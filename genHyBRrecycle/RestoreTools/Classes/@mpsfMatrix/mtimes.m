function y = mtimes(A, x)
%
%  Overload * operation for mpsfMatrix object.
%
%  Implements A * x and A' * x  for mpsfMatrix object
%  A and vector x.  Result is returned as a vector y. 
%

%  J. Nagy  2/22/12

if ( isa(A, 'mpsfMatrix') )

  boundary = A.boundary;
  switch boundary
    case 'zero'
      boundary = 0;
    case {'reflexive', 'reflective', 'neumann'}
      boundary = 'symmetric';
    case 'periodic'
      boundary = 'circular';
    otherwise
      error('Not a valid boundary condition')
  end

  M = A.matdata;
  [nrow_pad, ncol_pad, nframes_pad] = size(M);
  [nrowx, ncolx, nframesx] = size(x);

  padsize = [nrow_pad - nrowx, ncol_pad - ncolx]/2;
  if padsize ~= fix(padsize)
    error('Use only even sized images')
  end

  if A.transpose % Multiply A'*x
    
    xpad = padarray(x, padsize, boundary);
    ypad = real(ifft2(conj(M).*fft2(xpad)));
    ypad = ypad(padsize(1)+1:padsize(1)+nrowx, padsize(2)+1:padsize(2)+ncolx, :);
    y = sum(ypad, 3);

  else % Multiply A*x

    xpad = padarray(x, padsize, boundary);
    xpad = xpad(:,:,ones(1,nframes_pad));
    ypad = real(ifft2(M .* fft2(xpad)));    
    y = ypad(padsize(1)+1:padsize(1)+nrowx, padsize(2)+1:padsize(2)+ncolx, :);
  end

else

  error('incorrect argument type')

end
