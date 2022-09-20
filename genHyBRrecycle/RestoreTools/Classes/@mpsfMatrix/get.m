function val = get(A, prop_name)
%
%  GET Get mpsfMatrix properties from specified object and return
%      the value.
%
%     val = get(A, prop_name);
%
%  Valid choices for prop_name are:
%    'psf', 'matdata', 'type', 'boundary', 'transpose', 'center'
%

%  J. Nagy  2/16/12

switch prop_name
case 'psf'
  val = A.psf;
case 'matdata'
  val = A.matdata;
case 'type'
  val = A.type;
case 'boundary'
  val = A.boundary;
case 'transpose'
  val = A.transpose;
case 'center'
  val = A.center;
otherwise
  error([prop_name, 'Is not valid mpsfMatrix property'])
end
