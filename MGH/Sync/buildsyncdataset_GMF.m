function [sz] = buildsyncdataset_GMF(patient, seizure, dataPath, onset, offset)

% [SZ] = BUILDSYNCDATASET(PATIENT, SEIZURE, DATAPATH, EDFFILE, NS5FILE,
% ECOGSTARTTIME, ECOGENDTIME, ECOGCH, LFPCH) builds a structure SZ
% containing synchronized ECoG and LFP data for a given seizure of a given
% PATIENT. It also saves into .mat file in the DATAPATH/PATIENT directory.
% Requires a PATIENT.m file describing seizure informations (see SZINFO).
% ONSET and OFFSET describe the amount of time to add before and after
% seizure, resp.
% 
% TODO: !!Works for Neuroport patients only at the moment!!
%xi
% Example
%   [sz] = buildsyncdataset('MG49', 'Seizure36', '/Users/Shared/Data', 20, 20);
% 
% See also: BUILDECOGDATASET, SYNCECOGLFP
% 
% Author: Louis-Emmanuel Martinet <louis.emmanuel.martinet@gmail.com>

fprintf(['Working on ' patient ' ' seizure ': \n']);

info = szinfo(dataPath, patient, seizure);
info
if isempty(info)
    error(['No seizure information for ' patient seizure '.']);
end
ecogCh = info.ECoG.Channels;
lfpCh = info.LFP.Channels;
ecogSyncCh = info.ECoG.SyncCh;
lfpSyncCh = info.LFP.SyncCh;

[ecogRef, ecogProp] = openEDF(info.ECoG.EdfFile, ecogSyncCh); 
ecogFs = ecogProp.SampleRate(1);
ecogSzOn = max([1, round((info.StartTime - onset) * ecogFs)]);
ecogSzOff = min([round((info.EndTime + offset) * ecogFs), size(ecogRef,1)]);
lfpProp = NSX_open(info.LFP.Ns5File);
% lfpRef = NSX_read(lfpProp, lfpSyncCh, 1, 0, Inf)';
lfpFs = 1/lfpProp.Period;
lfpMaxIdx = length(lfpRef);
fprintf('Reference elec. loaded\n');

% Load only the ecogCh ECoG channels (there are also EEG data)
[dECoG, ecogProp] = openEDF(info.ECoG.EdfFile, ecogCh);
dECoG = dECoG(ecogSzOn : ecogSzOff, :);
ecogRef = ecogRef(ecogSzOn : ecogSzOff);
% fclose(ecogProp.FILE.FID);
fprintf('ECoG loaded');

% Do LFP/ECoG syncing:
[ecogIdx, lfpIdx, ecogRealFs] = syncecoglfp_GMF(ecogRef, ecogFs, lfpRef, lfpFs);

% find beginning of lfp time by working backwards from first sync in ecog
tECoG = (0:ecogSzOff-ecogSzOn)/ecogFs - onset;
dt_onset = onset + tECoG(ecogIdx(1));
tmn = lfpIdx(1)/lfpFs - tECoG(ecogIdx(1)); % find when tECoG = 0 occurs in lfp time
lfpSzOn = max([1, round((tmn - onset) * lfpFs)]);
lfpSzOff = min([round((tmn+info.EndTime-info.StartTime + offset) * lfpFs), lfpMaxIdx]);
tLFP = (0:lfpSzOff-lfpSzOn)/lfpFs - onset;

% for each sync event, skip ahead ~500 bins in tLFP (~18 ms) to match tECoG
%   (31.656 sec per sync in LFP vs. 31.674 per sync in ECoG)

ecogIdx0 = [1; ecogIdx];
lfpIdx0 = [1; lfpIdx - lfpSzOn];
count = 1;
lfp_t_ind = [];

for i = 1 : length(ecogIdx0)-1
    dt_ecog = tECoG(ecogIdx0(i+1))-tECoG(ecogIdx0(i));
    dt_lfp = tLFP(lfpIdx0(i+1))-tLFP(lfpIdx0(i));
    
    % compute time difference (in # of bins)
    dt_diff = dt_ecog - dt_lfp;
    dt_diff_bins = round(dt_diff * lfpFs);
    
    % clip given # of bins from LFP time axis
    dt_i = round(dt_lfp*lfpFs);
    lfp_t_ind = [lfp_t_ind, count : count + dt_i - 1];
    count = count + dt_i + dt_diff_bins;
end
lfp_t_ind = [lfp_t_ind count:lfpSzOff-lfpSzOn];

% trim time axis
tLFP = tLFP(lfp_t_ind);
fprintf('Synchronized indexing found\n');

% Load EEG (if it exists)
if isfield(info, 'EEG')
    if isfield(info.EEG, 'Labels')
        eegCh = label2index(info.EEG.Labels, ecogProp.Label);
        [dEEG, ~] = openEDF(info.ECoG.EdfFile, eegCh);
        dEEG = dEEG(ecogSzOn : ecogSzOff, :);
        fprintf('EEG loaded\n');
    end
end

% Now get the original LFP data
lfpProp = NSX_open(info.LFP.Ns5File);
dLFP = NSX_read(lfpProp, lfpCh(1), lfpCh(end), lfpSzOn, lfpSzOff - lfpSzOn + 1, 'p', 'p')';
%fclose(lfpProp.FID);
fprintf('LFP loaded');

% clip LFP to window of interest
dLFP = dLFP(1:length(lfp_t_ind),:);
% dLFP = dLFP(:);
fprintf(', clipped (synced)\n');

% Create structures with synced data and its meta-data
ecog = struct('Name', 'ecog', 'File', ecogProp.FileName, 'SamplingRate', ecogFs, ...
    'SamplingRate_Real', ecogRealFs, 'NbChannels', length(ecogCh), ...
    'StartTime', ecogSzOn / ecogFs, 'EndTime', ecogSzOff / ecogFs, ...
    'Labels', ecogProp.Label(ecogCh, :), 'PhysLim',  struct('PhysMin', ...
    ecogProp.PhysMin(ecogCh), 'PhysMax', ecogProp.PhysMax(ecogCh)), ...
    'Data', dECoG, 'Time', tECoG);
lfpLabels = cell2mat(arrayfun(@(i) num2str(i, '%02d'), lfpCh, 'UniformOutput', false)');
lfp = struct('Name', 'lfp', 'File', [lfpProp.Path lfpProp.Filename], 'SamplingRate', lfpFs, ...
    'NbChannels', length(lfpCh), 'StartTime', lfpSzOn / lfpFs, ...
    'EndTime', lfpSzOff / lfpFs, 'Labels', lfpLabels, ...
    'Data', dLFP, 'Time', tLFP);
sz = struct('Patient', patient, 'Seizure', seizure, 'Onset', onset, ...
            'Offset', offset, 'ECoG', ecog, 'LFP', lfp);

if isfield(info, 'EEG')
    if isfield(info.EEG, 'Labels')
        eeg = struct('Name', 'eeg', 'File', ecogProp.FileName, 'SamplingRate', ecogRealFs, ...
            'NbChannels', length(eegCh), 'StartTime', ecogSzOn / ecogFs, ...
            'EndTime', ecogSzOff / ecogFs, 'Labels', ecogProp.Label(eegCh, :), ...
            'PhysLim', struct('PhysMin', ecogProp.PhysMin(eegCh), 'PhysMax', ecogProp.PhysMax(eegCh)), ...
            'Data', dEEG);
        sz.EEG = eeg;
    end
end

% Create structure with syncing information
syncCheck = struct('ecogRef', ecogRef, 'lfpRef', lfpRef, 'ecogIdx', ecogIdx, 'lfpIdx', lfpIdx, 'lfpSzOn', lfpSzOn, 'lfpSzOff', lfpSzOff);
sz.sync = syncCheck;

% Save data
fprintf('Saving...');
% save([dataPath '/Synced_Multiscale_Mats/' patient '_' seizure '_LFP_ECoG_EEG.mat'], '-v7.3', 'sz');
save(['~/' patient '_' seizure '_LFP_ECoG_EEG.mat'], '-v7.3', 'sz');
fprintf('done.\n');

% Plot the syncing results.
% plot(tLFP,lfpRef/max(lfpRef),'b',tECoG,ecogRef,'r');

end
