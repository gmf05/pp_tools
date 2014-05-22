function strarray = cell2str(cell_in,dim)
% function strcell = str2cell(strs,dim)
% %%%%%%%%%
% str2cell.m takes an array of strings (strs) and returns a cell
% full of strings.
% By default, operates along rows of strs; if dim == 2, operates along
% columns
% %%%%%%%%%
%
if nargin<2, dim = 1; end; % operation dimension -- default: rows
N = length(cell_in);
max_str_len = 0;
for n = 1:N, max_str_len = max(max_str_len, length(cell_in{n})); end
strarray = char(N,max_str_len);

switch dim    
    case 1 % operate along rows
        for n=1:N, strarray(n,1:length(cell_in{n})) = cell_in{n}; end;        
    case 2 % operate along columns
      for n=1:N, strarray(1:length(cell_in{n}),n) = cell_in{n}; end; 
    otherwise
        error('unrecognized value -- dim should be 1 (rows) or 2 (columns)');
        
end

end
