function val = IRget(options,name,default,flag)
%
%   VAL = IRget(options,name,default,flag)
%
%   IRget gets options parameters.
%   VAL = IRget(OPTIONS,'NAME') extracts the value of the named parameter
%   from IR options structure OPTIONS, returning an empty matrix if
%   the parameter value is not specified in OPTIONS.  It is sufficient to
%   type only the leading characters that uniquely identify the
%   parameter.  Case is ignored for parameter names.  [] is a valid OPTIONS
%   argument.
%
%   VAL = IRget(OPTIONS,'NAME',DEFAULT) extracts the named parameter as
%   above, but returns DEFAULT if the named parameter is not specified
%   in OPTIONS.  For example
%
%     param = IRget(options,'Rtol',1e-6);
%
%   returns Rtol = 1e-6 if the Rtol property is not specified in opts.
%
%   See also IRset and IR methods
%
%   J. Nagy August, 2011
%

% undocumented usage for fast access with no error checking
if (nargin == 4) && isequal(flag,'fast')
    val = IRgetfast(options,name,default);
    return
end

if nargin < 2
    error('Not enough input arguments.');
end
if nargin < 3
    default = [];
end

if ~isempty(options) && ~isa(options,'struct')
    error('First argument must be an options structure created with HyBRset.');
end

if isempty(options)
    val = default;
    return;
end
allfields = {'x0'; 'MaxIter';'x_true';'Rtol';'IterBar'};

Names = allfields;

name = deblank(name(:)'); % force this to be a row vector
j = find(strncmpi(name,Names,length(name)));
if isempty(j)               % if no matches
    error(['Unrecognized property name ''%s''.  ' ...
        'See HyBRset for possibilities.'], name);
elseif length(j) > 1            % if more than one match
    % Check for any exact matches (in case any names are subsets of others)
    k = find(strcmpi(name,Names));
    if length(k) == 1
        j = k;
    else
        msg = sprintf('Ambiguous property name ''%s'' ', name);
        msg = [msg '(' Names{j(1),:}];
        for k = j(2:length(j))'
            msg = [msg ', ' Names{k,:}];
        end
        msg = sprintf('%s).', msg);
        error(msg);
    end
end

if any(strcmp(Names,Names{j,:}))
    val = options.(Names{j,:});
    if isempty(val)
        val = default;
    end
else
    val = default;
end

%------------------------------------------------------------------
function value = IRgetfast(options,name,defaultopt)
%IRGETFAST- Get IR OPTIONS parameter with no error checking.
%   VAL = IRGETFAST(OPTIONS,FIELDNAME,DEFAULTOPTIONS) will get the
%   value of the FIELDNAME from OPTIONS with no error checking or
%   fieldname completion. If the value is [], it gets the value of the
%   FIELDNAME from DEFAULTOPTIONS, another OPTIONS structure which is
%   probably a subset of the options in OPTIONS.
%
if isempty(options)
     value = defaultopt.(name);
     return;
end
% We need to know if name is a valid field of options, but it is faster to use 
% a try-catch than to test if the field exists and if the field name is
% correct.
try
    value = options.(name);
catch
    value = [];
    lasterr('');  % clean up last error
end

if isempty(value)
    value = defaultopt.(name);
end


