function varargout = callObj(name ,obj, args)

resultsize = nargout;
if nargin == 2
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
varargout = recfind(varargout);

end