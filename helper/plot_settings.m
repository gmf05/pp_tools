global PLOT_COLOR DO_CONF_INT DO_MASK

% plot settings 
PLOT_COLOR = 'b';
DO_CONF_INT = true;
DO_MASK = true;

% setting some default options:

% figure docking
% set(0,'DefaultFigureWindowStyle','docked') %dock the figures	
% set(0,'DefaultFigureWindowStyle','normal') %undock the figures
% set(0,'DefaultFigureWindowStyle','modal') %undock the figures

% figure sizing
% set(0,'DefaultFigureUnits','normalized');
% set(0,'DefaultFigurePosition',[0 0 1 1]);

% font preferences in figures
fontname = 'Helvetica';
fontsize = 12;
set(0,'defaultaxesfontname',fontname);
set(0,'defaulttextfontname',fontname);
set(0,'defaultaxesfontsize',fontsize);
set(0,'defaulttextfontsize',fontsize);


