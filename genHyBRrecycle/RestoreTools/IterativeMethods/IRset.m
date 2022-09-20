function options = IRset(varargin)
%
%   OPTIONS = IRset(varargin)
%
%   Create/alter options structure for iterative image restoration codes.
%   OPTIONS = IRset('PARAM1',VALUE1,'PARAM2',VALUE2,...) 
%     creates an options structure in which the named parameters have
%     the specified values.  Any unspecified parameters are set to [] 
%     (parameters with value [] indicate to use the default value for that 
%     parameter when passed to an IR function). It is sufficient to type
%     only the leading characters that uniquely identify the parameter.  
%     Case is ignored for parameter names.
%     NOTE: For values that are strings, the complete string is required.
%
%   OPTIONS = IRset(OLDOPTS,'PARAM1',VALUE1,...) creates a copy of
%     OLDOPTS with the named parameters altered with the specified values.
%
%   OPTIONS = IRset(OLDOPTS,NEWOPTS) combines an existing options structure
%     OLDOPTS with a new options structure NEWOPTS.  Any parameters in NEWOPTS
%     with non-empty values overwrite the corresponding old parameters in
%     OLDOPTS.
%
%   IRset with no input arguments and no output arguments displays all
%     parameter names and their possible values, with defaults shown in {}
%
%   OPTIONS = IRset(with no input arguments) creates an options structure
%     where all the fields are set to [].
%
%   OPTIONS = IRset('IR') creates an options structure with all
%     the parameter names and default values relevant to an IR method. 
%     That is,
%           IRset('IR')
%   or
%           IRset(@IR)
%   returns an options structure containing all the parameter names and
%   default values relevant to an IR method.
%
% IRset PARAMETERS for MATLAB: ( default parameter in {} )
%       x0 - initial guess
%               [ array | {none} ]
%   MaxIter - Maximum number of iterations 
%                   [positive integer | {[]} ]
%   x_true - True solution : [ array | {off} ]
%                Returns error with respect to x_true at each iteration
%  Rtol - Relative residual norm, norm(b-A*x)/norm(b), stopping tolerance
%               [ non-negative scalar | {10^-6}]
%  NE_Rtol - Relative residual norm, norm(b-A*x)/norm(b), stopping tolerance
%               [ non-negative scalar | {10^-6}]
%    IterBar - turn off/on waitbar that shows iteration progress
%                 [ off | {on} ]
%  tau - Step length for gradient descent methods
%                 [ scalar | {0} ]
%        Here, tau = 0 ==> steepest descent
%              tau < 0 ==> fixed, default value, results in Richardson 
%                                 iteration
%              tau > 0 ==> fixed with given value, again result in 
%                          Richardson iteration
%  TikParam - Tikhonov regularization parameter.  If alpha is not
%          specified, then it is assumed that alpha = 0, and
%          regularization needs to be enforced by the stopping 
%          criteria. Except in the case of IRhybr, which attempts
%          to choose a regularization parameter through GCV.
%  TikFun - Can be "Identity" or "Laplacian".  If not specified,
%           and TikParam is also not specified, then it is assumed
%           that regularization is enforced by the stopping iteration.
%           If TikParam is specified, but TikFun is not, then the
%           default is "Identity"
%
%   Examples:
%     To create OPTIONS with the default options for an IR nethod
%       OPTIONS = IRset('IR');
%     To change the maximum iterations to 150 in OPTIONS
%       OPTIONS = IRset(OPTIONS,'MaxIter',150);
%
%   See also Iterative Methods for Image Restoration
%
%   J. Nagy August, 2011

% Print out possible values of properties. 
if (nargin == 0) && (nargout == 0)
  fprintf('            x0: [ array | {none} ]\n');
  fprintf('           MaxIter: [ positive integer  | {[]} ]\n');
  fprintf('            x_true: [ array | {off} ]\n');
  fprintf('              Rtol: [  non-negative scalar | {10^-6}  ]\n');
  fprintf('           NE_Rtol: [  non-negative scalar | {10^-6}  ]\n');
  fprintf('           IterBar: [ off | {on} ]\n')
  fprintf('        StepLength: [ scalar | {0} ]\n');
  fprintf('            TikFun: [ {none} | Identity | Laplacian ]\n');
  fprintf('          TikParam: [ non-negative scalar | {gcv} ]\n' );
  return;
end

% Create a struct of all the fields
allfields = {'x0'; 'MaxIter';'x_true';'Rtol';'NE_Rtol';'IterBar'; ...
  'StepLength';'TikFun';'TikParam'};
  
% create cell array
structinput = cell(2,length(allfields));
% fields go in first row
structinput(1,:) = allfields';
% []'s go in second row
structinput(2,:) = {[]};
% turn it into correctly ordered comma separated list and call struct
options = struct(structinput{:});

numberargs = nargin; % we might change this value, so assign it


% If we pass in a function name then return the defaults.
if (numberargs==1) && (ischar(varargin{1}) || isa(varargin{1},'function_handle') )
  if ischar(varargin{1})
    funcname = lower(varargin{1});
      if ~exist(funcname)
        error('Undefined function.  Please use IRset(@funcname), where funcname is an IR method')
      end
  elseif isa(varargin{1},'function_handle')
    funcname = func2str(varargin{1});
  end
  try 
    optionsfcn = feval(varargin{1},'defaults');
  catch
    error('IRset ONLY works with an IR method.  Please use IRset(@IR).')
  end
  % The defaults from an IR functions don't include all the fields
  % typically, so run the rest of IRset as if called with
  % IRset(options,optionsfcn) to get all the fields.
  varargin{1} = options;
  varargin{2} = optionsfcn;
  numberargs = 2;
end

Names = allfields;
m = size(Names,1);
names = lower(Names);

i = 1;
while i <= numberargs
  arg = varargin{i};
  if ischar(arg)                         % arg is an option name
    break;
  end
  if ~isempty(arg)                      % [] is a valid options argument
    if ~isa(arg,'struct')
      error(['Expected argument %d to be a string parameter name ' ...
          'or an options structure \n created with IRset.'], i);
    end
    for j = 1:m
      if any(strcmp(fieldnames(arg),Names{j,:}))
        val = arg.(Names{j,:});
      else
        val = [];
      end
      if ~isempty(val)
        if ischar(val)
          val = lower(deblank(val));
        end
        checkfield(Names{j,:},val)
        options.(Names{j,:}) = val;
      end
    end
  end
  i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(numberargs-i+1,2) ~= 0
  error('Arguments must occur in name-value pairs.');
end

expectval = 0;                          % start expecting a name, not a value
while i <= numberargs
  arg = varargin{i};

  if ~expectval
    if ~ischar(arg)
      error('Expected argument %d to be a string parameter name.', i);
    end

    lowArg = lower(arg);
    j = strmatch(lowArg,names);
    if isempty(j)                       % if no matches
      error( sprintf('Invalid parameter name ''%s'' ', arg));
    elseif length(j) > 1                % if more than one match
      % Check for any exact matches (in case any names are subsets of others)
      k = strmatch(lowArg,names,'exact');
      if length(k) == 1
        j = k;
      else
        error(sprintf('Ambiguous parameter name ''%s'' ', arg));
      end
    end
    expectval = 1;                      % we expect a value next

  else
    if ischar(arg)
      arg = lower(deblank(arg));
    end
    checkfield(Names{j,:},arg);
    options.(Names{j,:}) = arg;
    expectval = 0;
  end
  i = i + 1;
end

if expectval
  error( 'Expected value for parameter ''%s''.', arg);
end


%----SUBFUNCTION---------------------------------------------
function checkfield(field,value)
%CHECKFIELD Check validity of structure field contents.
%   CHECKFIELD('field',V) checks the contents of the specified
%   value V to be valid for the field 'field'. 
%

% empty matrix is always valid
if isempty(value)
  return
end

% See if 'field' is a valid field.
validfield = true;
switch field
  case {'MaxIter'} % real positive integer
    [validvalue, errmsg] = PosInteger(field,value);
  case {'x0'} % numeric array, none
    [validvalue, errmsg] = x0type(field,value);
  case {'IterBar'} % off, on
    [validvalue, errmsg] = onOffType(field,value);
  case {'x_true'} % numeric array, off
    [validvalue, errmsg] = x_truetype(field,value);
  case {'Rtol'}% real non-negative scalar
    [validvalue, errmsg] = nonNegscalar(field,value);
  case {'NE_Rtol'}% real non-negative scalar
    [validvalue, errmsg] = nonNegscalar(field,value);
  case {'StepLength'} % real scalar
    [validvalue, errmsg] = realScalar(field,value);
  case {'TikFun'} % character string
    [validvalue, errmsg] = tikFunType(field,value);
  case {'TikParam'} % non-negative scalar, gcv
    [validvalue, errmesg] = tikParamType(field,value);
  otherwise
    validfield = false;  
    validvalue = false;
    errmsg = sprintf('Unrecognized parameter name ''%s''.', field);
end

if validvalue 
    return;
else
  error(errmsg)
end

%------------------------------------------------------------------------

function [valid, errmsg] = PosInteger(field,value)
% Any positive real integer
valid =  isreal(value) && isscalar(value) && (value > 0) && value == floor(value) ;

if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a positive real integer.',field);
else
  errmsg = '';
end

%-----------------------------------------------------------------------

function [valid, errmsg] = onOffType(field,value)
% One of these strings: on, off
valid =  ischar(value) && any(strcmp(value,{'on';'off'}));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''off'' or ''on''.',field);
else
  errmsg = '';
end

%-----------------------------------------------------------------------

function [valid, errmsg] = x_truetype(field,value)
% Either a numeric array or off
valid =  isnumeric(value) | (ischar(value) && any(strcmp(value,{'off'})));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a numeric array or ''off''.',field);
else
  errmsg = '';
end

%-----------------------------------------------------------------------

function [valid, errmsg] = x0type(field,value)
% Either a numeric array or none
valid =  isnumeric(value) | (ischar(value) && any(strcmp(value,{'none'})));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a numeric array or ''none''.',field);
else
  errmsg = '';
end


%------------------------------------------------------------------------

function [valid, errmsg] = nonNegscalar(field,value)
% Any real non-negative scalar
valid =  isreal(value) && isscalar(value) && (value >= 0);

if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-negative scalar.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

function [valid, errmsg] = realScalar(field,value)
% Any real  scalar
valid =  isreal(value) && isscalar(value);

if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-negative scalar.',field);
else
  errmsg = '';
end

function [valid, errmsg] = tikParamType(field,value)
% Either a numeric array or none
valid =  isreal(value) && isscalar(value) && (value >= 0) || (ischar(value) && any(strcmp(value,{'gcv'})));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a numeric array or ''none''.',field);
else
  errmsg = '';
end

function [valid, errmsg] = tikFunType(field,value)
% One of these strings: none, Identity, Laplacian
valid =  ischar(value) && any(strcmp(value,{'none';'identity';'laplacian'}));
if ~valid
  errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''none'' or ''Identity'' or ''Laplacian''.',field);
else
  errmsg = '';
end