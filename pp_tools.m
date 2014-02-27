clear;
fprintf('\n  Point Process Toolbox: v 1.0\n');
fprintf('Grant Fiddyment  Boston U, 01/2014\n \n');

% global variables:
global PP_TOOLS
global PLOT_COLOR
global FONT_SIZE
global DO_SAVE
global RESCALE_FUNC

PLOT_COLOR = 'b';
FONT_SIZE = 18;
DO_SAVE = true;
RESCALE_FUNC = 'exp';

% MGH-specific global variables
global DATA_DIR
global SPIKE_THRESH_ECOG
global SPIKE_THRESH_LFP
global SPIKE_THRESH_MUA
global MIN_REFRACT

SPIKE_THRESH_ECOG = 1;
SPIKE_THRESH_LFP = 1;
SPIKE_THRESH_MUA = NaN; % will need to test this
MIN_REFRACT = 1;

% add paths:
run('~/Code/matlab/get_gmf_env'); 
addpath(genpath(PP_TOOLS));
addpath(genpath(DATA_DIR));


