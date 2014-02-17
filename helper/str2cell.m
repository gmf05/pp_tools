function strcell = str2cell(strs,dim)
% function strcell = str2cell(strs,dim)
% %%%%%%%%%
% str2cell.m takes an array of strings (strs) and returns a cell
% full of strings.
% By default, operates along rows of strs; if dim == 2, operates along
% columns
% %%%%%%%%%
%
if nargin<2, dim = 1; end; % operation dimension -- default: rows
N = size(strs,dim);
strcell = cell(1,N);
f = @strtrim;
switch dim    
    case 1 % operate along rows
        for n=1:N, strcell{n} = f(strs(n,:)); end;
        
    case 2 % operate along columns
        for n=1:N, strcell{n} = f(strs(:,n)); end;        
        
    otherwise
        error('unrecognized value -- dim should be 1 (rows) or 2 (columns)');
        
end

end