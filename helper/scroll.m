function scroll(stepSize,numScrolls)
% scroll(stepSize,numScrolls) : Scrolls x-axis on a plot

% subplots
subplts = get(gcf,'Children');
Nsubplts = length(subplts);
for n = 1:numScrolls
  for s = 1:Nsubplts
    subplot(subplts(s))
    xlim(get(gca,'XLim')+stepSize);
  end
  if numScrolls>1, pause; end
end

% % % no subplots
% % for n = 1:numScrolls, xlim(get(gca,'XLim')+stepSize);  pause; end

end