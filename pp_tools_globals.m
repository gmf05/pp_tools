global PP_TOOLS RESCALE_FUNC PLOT_COLOR FONT_SIZE DO_SAVE DO_CONF_INT ...
  DO_MASK DATA

% path variables
PP_TOOLS = getenv('PP_TOOLS');
DATA = getenv('DATA');

% function used for rescaled interspike intervals
RESCALE_FUNC = 'exp';

% plot settings 
PLOT_COLOR = 'b';
FONT_SIZE = 18;
DO_SAVE = true;
DO_CONF_INT = true;
DO_MASK = true;


