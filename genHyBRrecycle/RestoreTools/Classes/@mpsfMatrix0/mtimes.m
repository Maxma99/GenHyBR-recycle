function y = mtimes(A, x)
%
%  Overload * operation for mpsfMatrix0 object.
%
%  Implements A * x and A' * x  for mpsfMatrix0 object
%  A and vector x.  Result is returned as a vector y. 
%

%  J. Nagy  2/23/12

if ( isa(A, 'mpsfMatrix0') )

  M = A.matdata;
  [m, n, nframes] = size(A.psf);
  xvec = (numel(x) == size(x,1));
  if A.transpose % Multiply A'*x
    if xvec
      x = reshape(x, m, n, nframes);
    end
    y = real(ifft2(conj(M).*fft2(x)));
    y = sum(y, 3);
    if xvec
      y = y(:);
    end
  else % Multiply A*x
    if xvec
      x = reshape(x, m, n);
    end
    x = x(:,:,ones(1,nframes));
    y = real(ifft2(M .* fft2(x)));
    if xvec
      y = y(:);
    end    
  end

else

  error('incorrect argument type')

end
