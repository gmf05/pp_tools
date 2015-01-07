% clear;
fprintf('\n  Point Process Toolbox: v 1.0\n');
fprintf('Grant Fiddyment  Boston U, 01/2014\n \n');

% add paths to all functions:
PP_TOOLS='~/Code/repos/pp_tools';
addpath(PP_TOOLS);
addpath([PP_TOOLS '/classes']);
addpath([PP_TOOLS '/helper']);
addpath([PP_TOOLS '/examples']);
addpath([PP_TOOLS '/python']);

 % initialize some global variables used for plotting:
plot_settings();

% function used for rescaled interspike intervals (can we get rid of this?)
% global RESCALE_FUNC
% RESCALE_FUNC = 'exp';