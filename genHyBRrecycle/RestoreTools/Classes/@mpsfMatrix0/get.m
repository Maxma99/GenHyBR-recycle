function val = get(A, prop_name)
%
%  GET Get mpsfMatrix0 properties from specified object and return
%      the value.
%
%     val = get(A, prop_name);
%
%  Valid choices for prop_name are:
%    'psf', 'matdata', 'transpose', 'center'
%

%  J. Nagy  2/23/12

switch prop_name
case 'psf'
  val = A.psf;
case 'matdata'
  val = A.matdata;
case 'transpose'
  val = A.transpose;
case 'center'
  val = A.center;
otherwise
  error([prop_name, 'Is not valid mpsfMatrix0 property'])
end
