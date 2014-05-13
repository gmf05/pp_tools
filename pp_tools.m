clear;
fprintf('\n  Point Process Toolbox: v 1.0\n');
fprintf('Grant Fiddyment  Boston U, 01/2014\n \n');

% global variables:
global PP_TOOLS
global PLOT_COLOR
global FONT_SIZE
global DO_SAVE
global RESCALE_FUNC
global DO_CONF_INT
global DO_MASK

PLOT_COLOR = 'b';
FONT_SIZE = 18;
DO_SAVE = true;
RESCALE_FUNC = 'exp';
DO_CONF_INT = true;
DO_MASK = true;

% MGH-specific global variables
global DATA_DIR
global SPIKE_THRESH_ECOG
global SPIKE_THRESH_LFP
global SPIKE_THRESH_MUA
global MIN_REFRACT

SPIKE_THRESH_ECOG = 0.8;
% SPIKE_THRESH_LFP = 0.4;
SPIKE_THRESH_LFP = 3;
SPIKE_THRESH_MUA = NaN; % will need to test this
MIN_REFRACT = 1;

% add paths:
PP_TOOLS = pwd();
% DATA_DIR = '/projectnb/ecog/Data'; % SCC
DATA_DIR = '/media/Shared/GMF/Documents/BostonU/Research/Data'; % SCC
addpath(genpath(PP_TOOLS));
addpath(genpath(DATA_DIR));


