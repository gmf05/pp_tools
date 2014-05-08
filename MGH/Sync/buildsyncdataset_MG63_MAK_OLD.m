function [sz] = buildsyncdataset_MG63_MAK(patient, seizure, dataPath, onset, offset)

% [SZ] = BUILDSYNCDATASET(PATIENT, SEIZURE, DATAPATH, EDFFILE, NS5FILE,
% ECOGSTARTTIME, ECOGENDTIME, ECOGCH, LFPCH) builds a structure SZ
% containing synchronized ECoG and LFP data for a given seizure of a given
% PATIENT. It also saves into .mat file in the DATAPATH/PATIENT directory.
% Requires a PATIENT.m file describing seizure informations (see SZINFO).
% ONSET and OFFSET describe the amount of time to add before and after
% seizure, resp.
% 
% TODO: !!Works for Neuroport patients only at the moment!!
%
% Example
%   [sz] = buildsyncdataset('MG49', 'Seizure36', '/Users/Shared/Data', 20, 20);
% 
% See also: BUILDECOGDATASET, SYNCECOGLFP
% 
% Author: Louis-Emmanuel Martinet <louis.emmanuel.martinet@gmail.com>

fprintf(['Working on ' patient ' ' seizure ': ']);

info = szinfo_MAK(dataPath, patient, seizure);
if isempty(info)
    error(['No seizure information for ' patient seizure '.']);
end
ecogCh = info.ECoG.Channels;
lfpCh = info.LFP.Channels;

ecogSyncCh = 129; % channel 129 in ECoG carries sync events
lfpSyncCh = 97; % channel 97 in Neuroport LFP carries sync events
[ecogRef, ecogProp] = openEDF(info.ECoG.EdfFile, ecogSyncCh); 
ecogFs = ecogProp.SampleRate(1);
%lfpProp = NSX_open(info.LFP.Ns5File);
%lfpRef = NSX_read(lfpProp, lfpSyncCh, 1, 0, Inf)';
%lfpFs = 1/lfpProp.Period;
d = openNSx('read', info.LFP.Ns5File, ['c:' num2str(lfpSyncCh) ':' num2str(lfpSyncCh)]);
lfpRef = d.Data';
lfpFs = d.MetaTags.SamplingFreq;
[ecogIdx, lfpIdx, ecogRealFs] = syncecoglfp(ecogRef, ecogFs, lfpRef, lfpFs);
lfpMaxIdx = length(lfpRef);
%clear ecogRef lfpRef;
fprintf('Synchronized indexes built');

% Load only the ecogCh ECoG channels (there are also EEG data) 
[dECoG, ecogProp] = openEDF(info.ECoG.EdfFile, ecogCh);
ecogSzOn = max([1, round((info.StartTime - onset) * ecogFs)]);
ecogSzOff = min([round((info.EndTime + offset) * ecogFs), size(dECoG,1)]);
dECoG = dECoG(ecogSzOn : ecogSzOff, :);
fclose(ecogProp.FILE.FID);
fprintf(', ECoG loaded');

% Find the sync events during the sz
idxSz = find(ecogSzOn <= ecogIdx & ecogIdx <= ecogSzOff);

% Time difference between ecogSzOn and the first sync event
dtStartSz = (ecogIdx(idxSz(1)) - ecogSzOn) / ecogRealFs;

% Use that to find lfpSzOn
lfpSzOn = round(lfpIdx(idxSz(1)) - dtStartSz * lfpFs);

% Same for the end
dtEndSz = (ecogSzOff - ecogIdx(idxSz(end))) / ecogRealFs;  %Last sync preceeding end of sz.
lfpSzOff = round(lfpIdx(idxSz(end)) + dtEndSz * lfpFs);

if lfpSzOn < 1 || lfpSzOff > lfpMaxIdx
    error('LFP data are truncated compared to ECoG data.');
end

% Now get the correct LFP data
%dLFP = NSX_read(lfpProp, lfpCh(1), lfpCh(end), lfpSzOn, lfpSzOff - lfpSzOn + 1, 'p', 'p')';
%fclose(lfpProp.FID);
dLFP = openNSx('read', info.LFP.Ns5File, ['c:' num2str(lfpCh(1)) ':' num2str(lfpCh(end))], ...
                                         ['t:' num2str(lfpSzOn)  ':' num2str(lfpSzOff)]);
dLFP = dLFP.Data';
fprintf(', LFP loaded');

ecog = struct('Name', 'ecog', 'File', ecogProp.FileName, 'SamplingRate', ecogRealFs, ...
    'NbChannels', length(ecogCh), 'StartTime', ecogSzOn / ecogFs, ...
    'EndTime', ecogSzOff / ecogFs, 'Labels', ecogProp.Label(ecogCh, :), ...
    'PhysLim', struct('PhysMin', ecogProp.PhysMin(ecogCh), 'PhysMax', ecogProp.PhysMax(ecogCh)), ...
    'Data', dECoG);
lfpLabels = cell2mat(arrayfun(@(i) num2str(i, '%02d'), lfpCh, 'UniformOutput', false)');
lfp = struct('Name', 'lfp', 'File', [info.LFP.Ns5File], 'SamplingRate', lfpFs, ...
    'NbChannels', length(lfpCh), 'StartTime', lfpSzOn / lfpFs, ...
    'EndTime', lfpSzOff / lfpFs, 'Labels', lfpLabels, ...
    'Data', dLFP);
sz = struct('Patient', patient, 'Seizure', seizure, 'Onset', onset, ...
            'Offset', offset, 'ECoG', ecog, 'LFP', lfp);

if isfield(info.EEG, 'Labels')
    eegCh = label2index(info.EEG.Labels, ecogProp.Label);
    [dEEG, ~] = openEDF(info.ECoG.EdfFile, eegCh);
    dEEG = dEEG(ecogSzOn : ecogSzOff, :);
    fprintf(', EEG loaded');
    eeg = struct('Name', 'eeg', 'File', ecogProp.FileName, 'SamplingRate', ecogRealFs, ...
    'NbChannels', length(eegCh), 'StartTime', ecogSzOn / ecogFs, ...
    'EndTime', ecogSzOff / ecogFs, 'Labels', ecogProp.Label(eegCh, :), ...
    'PhysLim', struct('PhysMin', ecogProp.PhysMin(eegCh), 'PhysMax', ecogProp.PhysMax(eegCh)), ...
    'Data', dEEG);
    sz.EEG = eeg;
end

syncCheck = struct('ecogRef', ecogRef, 'lfpRef', lfpRef, 'ecogIdx', ecogIdx, 'lfpIdx', lfpIdx);
sz.sync = syncCheck;

fprintf(', saving');
save([dataPath '/Synced_Multiscale_Mats/' patient '_' seizure '_LFP_ECoG_EEG.mat'], '-v7.3', 'sz');
fprintf(', done.\n');

% % Quick check plot
% edt = 1 / sz.ECoG.SamplingRate;
% e = sz.ECoG.Data;
% ldt = 1 / sz.LFP.SamplingRate;
% l = sz.LFP.Data;
% ax = plotyy(0:edt:(length(e)-1)*edt, e(:,35), 0:ldt:(length(l)-1)*ldt, l(:,48));
% legend('ECoG Channel 35', 'LFP Channel 48');
% xlabel('time (s)');
% set(get(ax(1),'Ylabel'),'String','Amplitude') 
% set(get(ax(2),'Ylabel'),'String','Amplitude') 

end

% load /Users/Shared/Data/MG49/MG49_Seizure45.mat
% eegFile = '/Users/Shared/Data/MG49/MG49_Notes_Images/MG49_eeg_channels.txt';
% fid = fopen(eegFile, 'r');
% eegLabels = textscan(fid, '%s'); eegLabels = eegLabels{1};
% [ecogRef, ecogProp] = openEDF('/Users/Shared/Data/MG49/MG49_edfs/MG49_Seizure36.edf', 1);
% cellLabels = cellstr(ecogProp.Label);
% eegCh = zeros(1, size(eegLabels, 1));
% for i = 1 : size(eegLabels, 1)
% func = @(s) strcmp(s, eegLabels{i});
% eegCh(i) = find(cellfun(func, cellLabels), 1);
% end
% ecogName = '/Users/Shared/Data/MG49/MG49_edfs/MG49_Seizure36.edf';
% [dEEG, ~] = openEDF(ecogName, eegCh);
% ecogSzOn = sz.ECoG.StartTime * sz.ECoG.SamplingRate;
% ecogSzOff = sz.ECoG.EndTime * sz.ECoG.SamplingRate;
% dEEG = dEEG(ecogSzOn : ecogSzOff, :);
% fprintf('EEG loaded\n');
% ecogRealFs = sz.ECoG.SamplingRate;
% eeg = struct('Name', 'eeg', 'File', ecogProp.FileName, 'SamplingRate', ecogRealFs, ...
% 'NbChannels', length(eegCh), 'StartTime', ecogSzOn / ecogRealFs, ...
% 'EndTime', ecogSzOff / ecogRealFs, 'Labels', ecogProp.Label(eegCh, :), ...
% 'Data', dEEG);
% sz.EEG = eeg;
% save /Users/Shared/Data/MG49/MG49_Seizure45.mat -v7.3 sz
