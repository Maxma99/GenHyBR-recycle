function A = set(A, varargin)
%
% Set mpsfMatrix0 properties and return the updated object.
%
%     A = set(A, ...)
% 
%  This funtion accepts an mpsfMatrix0 object, A, and a variable list of
%  property name/value pairs and returns the modified object.
%  Valid properties are:
%    'psf', 'matdata', 'transpose', 'center'
%

%  J. Nagy  2/23/12

property_argin = varargin;

while length( property_argin ) >= 2
  prop = property_argin{1};
  val = property_argin{2};
  property_argin = property_argin(3:end);

  switch prop
  case 'psf'
    A.psf = val;
  case 'matdata'
    A.matdata = val;
  case 'transpose'
    A.transpose = val;
  case 'center'
    A.center = val;
  otherwise
    error('Valid mpsfMatrix0 properties: psf, matdata, type, boundary, transpose')
  end
end

