function dCell = szinfo(sharedPath, patient, seizure)

% DCELL = SZINFO(SHAREDPATH, PATIENT[, SEIZURE]) loads seizure information
% if available. Returns the list of seizures and their infos if not
% particular SEIZURE are required. Requires a PATIENT.m file describing
% seizure informations. The seizure timing is defined relatively to an
% edfFile. Some constants are defined in this function and can be used (or
% redefined) in the PATIENT.m file: 
% - EDFPATH = sharedPath/patient/patient_edfs
% - LFPPATH = sharedPath/patient/patient_Neuroport
% - MAPPATH = sharedPath/patient/patient_Notes_Images
% 
% Example
%   info = szinfo('/Users/Shared/Data/', 'MG29', 'Seizure13');
% 
% Author: Louis-Emmanuel Martinet <louis.emmanuel.martinet@gmail.com>

patientDataFile = fullfile(sharedPath, patient, [patient '.m']);
dCell = [];
if ~exist(patientDataFile, 'file')
    return
end

EDFPATH = fullfile(sharedPath, patient, [patient '_edfs']);
MAPPATH = fullfile(sharedPath, patient, [patient '_Notes_Images']);
LFPPATH = fullfile(sharedPath, patient, [patient '_Neuroport']); %#ok<NASGU>

% get_gmf_env; SZSHAREDPATH = DATA_DIR;
% ecogSyncCh = 1; lfpSyncCh = 97;
eval(['run ''' patientDataFile '''']);
if nargin == 3
    szNames = {dCell.Seizure};
    func = @(s) strcmp(s, seizure);
    dCell = dCell(cellfun(func, szNames));
end
for i = 1 : length(dCell)
    if ~isfield(dCell(i).ECoG, 'EdfFile')
        dCell(i).ECoG.EdfFile = fullfile(EDFPATH, [patient '_' dCell(i).Seizure '.edf']); %#ok<AGROW>
    end
    if ~isfield(dCell(i).ECoG, 'MapFile')       
        dCell(i).ECoG.MapFile = fullfile(MAPPATH, [patient '_ecog.xls']); %#ok<AGROW>
    end
    if ~isfield(dCell(i), 'MatFile') || isempty(dCell(i).MatFile)
        dCell(i).MatFile = fullfile(sharedPath, patient, [patient '_' dCell(i).Seizure '.mat']); %#ok<AGROW>
    end    
    if isfield(dCell(i),'LFP') && isfield(dCell(i).LFP, 'Ns5File') && ~isfield(dCell(i).LFP, 'MapFile')
        dCell(i).LFP.MapFile = fullfile(MAPPATH, [patient '_lfp.xls']); %#ok<AGROW>
    end
    if isfield(dCell(i),'ECoG') && exist('ecogSyncCh','var')
        dCell.ECoG.SyncCh = ecogSyncCh;
    end
    if isfield(dCell(i),'LFP') && exist('lfpSyncCh','var')
        dCell.LFP.SyncCh = lfpSyncCh;
    end
end

end
