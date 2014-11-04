function print_latex_table(table,title)
% print_latex_table.m 
% Takes a MATLAB cell as input (table) and writes LaTeX syntax around
% the entries for pretty printing
%

if nargin<2, title = ''; end

% get size of table
N_row = size(table,1);
N_col = size(table,2);


% alignments:
table_alignment = 'center'; % (center, left, right)
col_alignment = 'c'; % (c, l, r)

% lines separating data:
row_line = '\\hline';
col_line = '|';

% % % use these settings to eliminate row lines,
% % % column lines, or both
row_line = '';
col_line = ''; 

% clc;
fprintf('\n');

% command to initialize table
fprintf(['\\begin{' table_alignment '} ' title '\\\\\n\\begin{tabular}{|']);
for j=1:N_col, fprintf([col_alignment col_line]); end;
if isempty(col_line)
    fprintf(['|} \\hline \n']);
else
    fprintf(['} \\hline \n']);
end

% useful for later:
end_row_text = ['\\\\ ' row_line ' \n'];

for i=1:N_row
    
    for j=1:N_col
                
        if ischar(table{i,j})
            str = table{i,j};
            str = strrep(str,'\','\\'); % also, replace \ by \\
        elseif isfloat(table{i,j})
            str = num2str(table{i,j},3); % NOTE: second argument governs printing
        else
            % any other cases?
            error('unknown data type');
        end
        fprintf(str);
        
        if j<N_col, fprintf(' & '); end
        
    end
    
    if i<N_row, fprintf(end_row_text); end
end


fprintf(['\\\\ \\hline \n \\end{tabular} \\end{' table_alignment '} \n \n \n']);

end

% reset float precision?


