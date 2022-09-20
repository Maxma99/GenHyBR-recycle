function varargout = size( A, dim )
%
%  Overload size for mpsfMatrix0 object.
%  This just displays the storage requirements for the matrix
%  data.
%

% J. Nagy  2/23/12

d1 = size(A.psf);
if length(d1) == 2
  d1(3) = 1;
end
d = [d1(3)*d1(1)^2, d1(2)^2];

if nargin == 2
  if dim <= length(d)
    d = d(dim);
  else
    d = 1;
  end
end

if nargout == 1 || nargout == 0
  varargout{1} = d;
else
  for i = 1:length(d)
    varargout{i} = d(i);
  end
end
