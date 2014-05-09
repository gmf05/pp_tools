function indexes = label2index(labels, labelList, noError)

% INDEXES = LABEL2INDEX(LABELS, LABELLIST, NOERROR) converts a set of
% LABELS to their indexes based on the full list of labels LABELLIST. If
% NOERROR = true, no error message is displayed when a label is not found.
% 
% Examples:
%   sz = loadsz('MG49', 'Seizure36');
%   id = label2index({'LATD1','SbFS4'}, sz.ECoG.Labels);

if nargin == 2
    noError = false;
end

if isempty(labels)
    indexes = [];
    return
end

if ~iscell(labelList)
    labelList = cellstr(labelList);
end

if ~iscell(labels)
    labels = cellstr(labels);
end

indexes = zeros(1, length(labels));
parfor i = 1 : length(labels)
    func = @(s) strcmp(s, labels{i});
    id = find(cellfun(func, labelList), 1);
    if isempty(id) && ~noError
        error(['Label ' labels{i} ' not found']);
    else
        indexes(i) = id;
    end
end

end