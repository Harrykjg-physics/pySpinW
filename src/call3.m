function varargout = call3(name ,obj, args, n)

if nargin < 4
    resultsize = nargout;
else
    resultsize = n;
end
if nargin < 3
    args = {};
end

if resultsize > 0
    % call the function with the given number of
    % output arguments:
    varargout = cell(resultsize, 1);
    if isempty(obj)
        [varargout{:}] = feval(name, args{:});
    else
        [varargout{:}] = feval(name, obj, args{:});
    end
else
    % try to get output from ans:
    clear('ans');
    if isempty(obj)
        feval(name, args{:})
    else
        feval(name, obj, args{:})
    end
    try
        varargout = {ans};
    catch err
        varargout = {[]};
    end
end

% Remove all non mapable reults
varargout = {varargout};
varargout = recfind(varargout);

end