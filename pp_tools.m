clear;
fprintf('\n  Point Process Toolbox: v 1.0\n');
fprintf('Grant Fiddyment  Boston U, 01/2014\n \n');

global PP_TOOLS RESCALE_FUNC PLOT_COLOR FONT_SIZE DO_SAVE DO_CONF_INT ...
  DO_MASK DATA

% add point process functions to MATLAB path
% add environment variable PP_TOOLS or hard-code path below...
PP_TOOLS = getenv('PP_TOOLS'); addpath(genpath(PP_TOOLS));
DATA = getenv('DATA'); addpath(DATA); % where MGH data is found

% function used for rescaled interspike intervals
RESCALE_FUNC = 'exp';

% plot settings 
PLOT_COLOR = 'b';
FONT_SIZE = 18;
DO_SAVE = true;
DO_CONF_INT = true;
DO_MASK = true;


