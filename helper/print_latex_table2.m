function print_latex_table2(table,title,filename)
% print_latex_table2.m
% NOTE TO SELF: CHANGE ALL PRINTING TO CONSOLE TO PRINTING TO FILE!!!
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

% printing to file instead of screen
texFile = [filename '.tex'];
delete(texFile);
fid = fopen(texFile,'w');

% preamble
fprintf(fid,'\\documentclass{article}\n');
% fprintf(fid,'\\usepackage{/home/gmf/Documents/gfiddy}\n');
fprintf(fid,'\\usepackage{amsmath}\n');
fprintf(fid,'\\setlength{\\parindent}{0pt}\n\n');
fprintf(fid,'\\begin{document}\n');

% command to initialize table
fprintf(fid,['\\begin{' table_alignment '} ' title '\\\\\n\\begin{tabular}{|']);
for j=1:N_col, fprintf(fid,[col_alignment col_line]); end;
if isempty(col_line)
    fprintf(fid,['|} \\hline \n']);
else
    fprintf(fid,['} \\hline \n']);
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
        fprintf(fid,str);
        
        if j<N_col, fprintf(fid,' & '); end
        
    end
    
    if i<N_row, fprintf(fid,end_row_text); end
end


fprintf(fid,['\\\\ \\hline \n \\end{tabular} \\end{' table_alignment '} \n \n \n']);

% end of document, more latex formatting
fprintf(fid,'\\end{document}\n');

% TeX -> PDF
eval(['!pdflatex ' texFile]);

end
