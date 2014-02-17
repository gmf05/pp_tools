function update_fig()
  
  global FONT_SIZE
  
  % get current figure
  my_fig = gcf(); 
    
  % for each subplot, change FONT_SIZE
  ax = findall(my_fig,'type','Axes');
  for i = 1:length(ax)
    lbl = get(ax(i),'XLabel');     
    set(lbl,'interpreter','latex');
    set(lbl,'FontSize',FONT_SIZE);
    
    lbl = get(ax(i),'YLabel');     
    set(lbl,'interpreter','latex');
    set(lbl,'FontSize',FONT_SIZE);
    
    lbl = get(ax(i),'Title');    
    set(lbl,'interpreter','latex');
    set(lbl,'FontSize',FONT_SIZE);
  end
  
end
