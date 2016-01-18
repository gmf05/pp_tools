% clear;
fprintf('\n  Point Process Toolbox: v 1.0\n');
fprintf('Grant Fiddyment  Boston U, 01/2014\n \n');

% add paths to all functions:
PP_TOOLS='~/Code/git/pp_tools';
addpath(PP_TOOLS);
addpath([PP_TOOLS '/classes']);
addpath([PP_TOOLS '/helper']);
addpath([PP_TOOLS '/examples']);
addpath([PP_TOOLS '/python']);

% function used for rescaled interspike intervals (can we get rid of this?)
% global RESCALE_FUNC
% RESCALE_FUNC = 'exp';

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

% line width
 set(0,'defaultlinelinewidth',3)
set(0,'defaultaxeslinewidth',2)
set(0,'defaultpatchlinewidth',2)

% font preferences in figures
fontname = 'Helvetica';
fontsize = 16;
set(0,'defaultaxesfontname',fontname);
set(0,'defaulttextfontname',fontname);
set(0,'defaultaxesfontsize',fontsize);
set(0,'defaulttextfontsize',fontsize);

