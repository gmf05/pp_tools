function SeizureGUI(varargin)

% declare some global variables
global DATA_DIR
global FIGURE
global MENUS
global BUTTONS
global PATIENT_LIST

% initialize figure
close all;
FIGURE = figure('units','normalized', ...
  'position', [0.01 0.01 0.9 0.9], ...
  'menubar', 'none', ...
  'name', 'SeizureGUI', ...
  'numbertitle', 'off', ...
  'resize', 'off');

% % % parse input??

% set path & patient list
DATA_DIR = '~/Documents/BostonU/Research/Data';
PATIENT_LIST = {'BU1'; 'BU3'; 'BU4'; 'BU7'; 'BU8'; ...
  'BU9'; 'BU12'; 'BU14'; 'BW09'; 'BW11'; 'BW17'; ...
  'BW19'; 'BW21'; 'MG08'; 'MG11'; 'MG13'; 'MG29'; 'MG31'; ...
  'MG32'; 'MG34'; 'MG37'; 'MG38b'; 'MG39b'; 'MG42'; ...
  'MG48'; 'MG49'; 'MG51'; 'MG56'; 'MG58'; 'MG61'; 'MG63'};

% create pull-down menus
MENUS(1) = uicontrol('style','pop',...
                 'unit','normalized',...
                 'position',[0.01 0.79 0.2 0.2],...
                 'backgroundc',get(FIGURE,'color'),...
                 'fontsize',12,'fontweight','bold',...
                 'string',PATIENT_LIST,'value',1,...
                 'callback',{@updatePatient});
MENUS(2) = uicontrol('style','pop',...
                 'unit','normalized',...
                 'position',[0.21 0.79 0.2 0.2],...
                 'backgroundc',get(FIGURE,'color'),...
                 'fontsize',12,'fontweight','bold',...
                 'callback',{@updateSeizure});
MENUS(3) = uicontrol('style','pop',...
                 'unit','normalized',...
                 'position',[0.41 0.79 0.2 0.2],...
                 'backgroundc',get(FIGURE,'color'),...
                 'fontsize',12,'fontweight','bold',...
                 'callback',{@updateData});
BUTTONS(1) = uicontrol('style','pushbutton',...
                 'unit','normalized',...
                 'position',[0.65 0.95 0.08 0.05],...
                 'backgroundc',get(FIGURE,'color'),...
                 'fontsize',12,'fontweight','bold',...
                 'string','Coord',...
                 'callback',{@getCoord});
BUTTONS(2) = uicontrol('style','pushbutton',...
                 'unit','normalized',...
                 'position',[0.73 0.95 0.05 0.05],...
                 'backgroundc',get(FIGURE,'color'),...
                 'fontsize',12,'fontweight','bold',...
                 'string','Edit',...
                 'callback',{@editPlotParams});
BUTTONS(3) = uicontrol('style','pushbutton',...
                 'unit','normalized',...
                 'position',[0.78 0.95 0.05 0.05],...
                 'backgroundc',get(FIGURE,'color'),...
                 'fontsize',12,'fontweight','bold',...
                 'string','Save',...
                 'callback',{@savePlot});
BUTTONS(4) = uicontrol('style','pushbutton',...
                 'unit','normalized',...
                 'position',[0.83 0.95 0.05 0.05],...
                 'backgroundc',get(FIGURE,'color'),...
                 'fontsize',12,'fontweight','bold',...
                 'string','Refresh',...
                 'callback',{@drawMontage});
               

% initialize seizure list & montage
updatePatient(MENUS(1),[]);

end

function updatePatient(src,evnt)

global PATIENT
global PATIENT_LIST

PATIENT = PATIENT_LIST{get(src,'Value')};
updateSeizureList();
updateDataList();
drawMontage();
fprintf(['\n\nPatient is ' PATIENT '\n\n']);

end

function updateSeizure(src,evnt)

global SEIZURE
global SEIZURE_LIST

SEIZURE = SEIZURE_LIST{get(src,'Value')};
updateDataList();
fprintf(['\n\nSeizure is ' SEIZURE '\n\n']);

end

function updateData(src,evnt)

global DATA_TYPE
global DATA_LIST

DATA_TYPE = DATA_LIST{get(src,'Value')};
fprintf(['\n\nData type is ' DATA_TYPE '\n\n']);

end

function updateSeizureList()

global PATIENT
global SEIZURE
global SEIZURE_LIST
global DATA_DIR
global MENUS

sz_list_full = char(fread(fopen([DATA_DIR '/SeizureList.txt'])))';
temp = sz_list_full;
[~,i] = regexp(temp,[PATIENT ':\n']);
temp = temp(i+1:end);
i = regexp(temp,'\n\n','once');
temp = temp(1:i);
nl_ind = regexp(temp,'\n');
N = length(nl_ind);
nl_ind = [0 nl_ind];

SEIZURE_LIST = cell(1,N);
for n = 1:N
  SEIZURE_LIST{n} = temp(nl_ind(n)+1:nl_ind(n+1)-1);
end

SEIZURE = SEIZURE_LIST{1};
set(MENUS(2),'string',SEIZURE_LIST,'value',1);

end

function updateDataList()

global PATIENT
global SEIZURE
global DATA_DIR
global DATA_LIST
global DATA_TYPE
global MENUS

if exist([DATA_DIR '/' PATIENT '/' PATIENT '_' SEIZURE '_LFP_ECoG_EEG.mat'],'file')
  DATA_LIST = {'ECoG'; 'LFP'};
else
  DATA_LIST = {'ECoG'};
end
DATA_TYPE = DATA_LIST{1};
set(MENUS(3),'string',DATA_LIST,'value',1);

end

function drawMontage(Ws)

cla; % clear axis

global PATIENT
global DATA_TYPE
global DATA_DIR
global CAXIS

if nargin<1 || range(Ws)==0
  Weights = false;
  col = 'r';
else
  Weights = true;
  colormap('default');
  if isempty(CAXIS), CAXIS = [min(Ws),max(Ws)]; end
  caxis(CAXIS);
  color_RGB = colormap();
end

% circle-drawing parameters:
R = 5; % radius
theta = 0:0.01:2*pi; % array of angles

switch DATA_TYPE
  case 'ECoG'
    PATIENT_DIR = [DATA_DIR '/' PATIENT];
    patient_img = imread([PATIENT_DIR '/' PATIENT '.png']);
    imagesc(patient_img), hold on;
    colormap(gray); freezeColors;
    load([PATIENT_DIR '/' PATIENT '_ECoG_map']);
    N_electrodes = length(ECoG.Coord);
  case 'LFP'
    [];
end

for n = 1:N_electrodes

  if Weights
    [~,W_ind] = min(abs(CAXIS - Ws(n)));
    col = color_RGB(W_ind,:);
  end

  x = ECoG.Coord(n,1);
  y = ECoG.Coord(n,2);
  xs = x+R*cos(theta');
  ys = y+R*sin(theta');
  fill(xs,ys,col);
end

end

function getCoord(src,evnt)

fig = figure('units','normalized', ...
  'position', [0.01 0.01 0.9 0.9], ...
  'menubar', 'none', ...
  'name', 'Get coordinates', ...
  'numbertitle', 'off', ...
  'resize', 'off');

end

function editPlotParams(src,evnt)

global CAXIS
% global WEIGHTS

fig = figure('units','normalized', ...
  'position', [0.01 0.4 0.4 0.5], ...
  'menubar', 'none', ...
  'name', 'Plot parameters', ...
  'numbertitle', 'off', ...
  'resize', 'off');
uicontrol('style','edit',...
           'unit','normalized',...
           'position',[0.01 0.79 0.2 0.2],...
           'backgroundc',get(fig,'color'),...
           'fontsize',12,'fontweight','bold',...
           'callback',{@updateWeights});
uicontrol('style','edit',...
           'unit','normalized',...
           'position',[0.01 0.5 0.2 0.2],...
           'backgroundc',get(fig,'color'),...
           'fontsize',12,'fontweight','bold',...
           'callback',{@updateCAxis});

% uicontrols for caxis, weights
% add button to close & drawMontage();

end

function updateWeights(src,evnt)
  eval(['Ws = ' get(src,'string') ';']);
  global FIGURE
  figure(FIGURE);
  drawMontage(Ws);
end

function updateCAxis(src,evnt)
  global CAXIS
  eval(['CAXIS = ' get(src,'string') ';']);
end

function savePlot(src,evnt)
end
