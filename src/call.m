function varargout = call(name, args)

resultsize = nargout;
if nargin == 1
    args = {};
end

if resultsize > 0
    % call the function with the given number of
    % output arguments:
    varargout = cell(resultsize, 1);
    [varargout{:}] = feval(name, args{:});
else
    % try to get output from ans:
    clear('ans');
    feval(name, args{:})
    try
        varargout = {ans};
    catch err
        varargout = {[]};
    end
end
end