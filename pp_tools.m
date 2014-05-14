clear;
fprintf('\n  Point Process Toolbox: v 1.0\n');
fprintf('Grant Fiddyment  Boston U, 01/2014\n \n');

global PP_TOOLS RESCALE_FUNC PLOT_COLOR FONT_SIZE DO_SAVE DO_CONF_INT ...
  DO_MASK

% add point process functions to MATLAB path
% add environment variable PP_TOOLS or hard-code path below...
PP_TOOLS = getenv('PP_TOOLS'); addpath(genpath(PP_TOOLS));

% function used for rescaled interspike intervals
RESCALE_FUNC = 'exp';

% plot settings 
PLOT_COLOR = 'b';
FONT_SIZE = 18;
DO_SAVE = true;
DO_CONF_INT = true;
DO_MASK = true;

% ============== MGH specific variables, paths ===================
% (1) where data is found -- hardcode if desired
global DATA; DATA = getenv('DATA'); addpath(DATA);

% (2) parameters for spike-finding
global SPIKE_THRESH_ECOG SPIKE_THRESH_LFP SPIKE_THRESH_MUA MIN_REFRACT
SPIKE_THRESH_ECOG = 0.8;
% SPIKE_THRESH_LFP = 0.4;
SPIKE_THRESH_LFP = 3;
SPIKE_THRESH_MUA = NaN; % will need to test this
MIN_REFRACT = 1;





