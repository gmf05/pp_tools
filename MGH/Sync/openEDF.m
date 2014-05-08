function [channelData, edfProperties] = openEDF(varargin)

% OPENEDF loads data from an edf file
% 
% [CHANNELDATA, EDFPROPERTIES] = OPENEDF([FILENAME,] [CH,] [TIMESTART,
% TIMEEND]) loads an edf file FILENAME. A subset of channels CH can be
% given, as well as a time interval [TIMESTART, TIMEEND] in sec to load. If
% no FILENAME is given, a window will open to select one. It results
% CHANNELDATA the matrix of data N x M (N time steps, M channels) and the
% properties EDFPROPERTIES of the recording session.
%
% Example
% [channelData, edfProperties] = openedf('MG49_Seizure36.edf', 2:21, 123, 149);
% 
% See also: SDFOPEN, SDFREAD
% 
% Author: Louis-Emmnuel Martinet (3/2012), louis.emmanuel.martinet@gmail.com

ch = [];
timeStart = []; 
timeEnd = [];
filePath = [];
if nargin == 1 && strcmp(class(varargin{1}), 'double')
    ch = varargin{1};
elseif nargin == 1
    filePath = varargin{1};
elseif nargin == 2 && strcmp(class(varargin{1}), 'double')
    timeStart = varargin{1};
    timeEnd = varargin{2};
elseif nargin == 2
    filePath = varargin{1};
    ch = varargin{2};
elseif nargin == 3 && strcmp(class(varargin{1}), 'double')
    ch = varargin{1};
    timeStart = varargin{2};
    timeEnd = varargin{3};    
elseif nargin == 3
    filePath = varargin{1};
    timeStart = varargin{2};
    timeEnd = varargin{3};
elseif nargin == 4
    filePath = varargin{1};
    ch = varargin{2};
    timeStart = varargin{3};
    timeEnd = varargin{4};
elseif nargin > 4
    fprintf('wrong number of arguments.\n');
    return;
end

if isempty(filePath)
    [filename, pathname, ~] = uigetfile('*.edf', 'Pick a Neural Recording');
    filePath = [pathname filename];
end

if isempty(ch)
    edfProperties = sdfopen(filePath, 'r');
else
    edfProperties = sdfopen(filePath, 'r', ch);
end

if isempty(timeStart)
    [channelData, edfProperties] = sdfread(edfProperties, Inf);
else
    [channelData, edfProperties] = sdfread(edfProperties, timeEnd - timeStart, timeStart);
end

end