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

fprintf(['Working on ' patient ' ' seizure '... \n']);

tic
info = szinfo(dataPath, patient, seizure);
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
ecogRef = ecogRef(ecogSzOn:ecogSzOff);

tic
lfpProp = NSX_open(info.LFP.Ns5File);
d = openNSx('read', info.LFP.Ns5File, ['c:' num2str(lfpSyncCh) ':' num2str(lfpSyncCh)], 'precision', 'double');
lfpRef = d.Data';
lfpFs = d.MetaTags.SamplingFreq;
lfpMaxIdx = length(lfpRef);
fprintf('Reference elec. loaded\n');
toc
  
% Load only the ecogCh ECoG channels (there are also EEG data)
tic
[dECoG, ecogProp] = openEDF(info.ECoG.EdfFile, ecogCh);
dECoG = dECoG(ecogSzOn : ecogSzOff, :);
fclose(ecogProp.FILE.FID);
fprintf('ECoG loaded\n');
toc

% Do LFP/ECoG syncing:
tic
[ecogIdx, lfpIdx, ecogRealFs] = syncecoglfp(ecogRef, ecogFs, lfpRef, lfpFs);

% find beginning of lfp time by working backwards from first sync in ecog
tECoG = (0:ecogSzOff-ecogSzOn)/ecogFs - onset;
dt_onset = onset + tECoG(ecogIdx(1));
tmn = lfpIdx(1)/lfpFs - tECoG(ecogIdx(1)); % find when tECoG = 0 occurs in lfp time
lfpSzOn = max([1, round((tmn - onset) * lfpFs)]);
lfpSzOff = min([round((tmn+info.EndTime-info.StartTime + offset) * lfpFs), lfpMaxIdx]);
tLFP = (0:lfpSzOff-lfpSzOn)/lfpFs - onset;
lfpIdx0 = lfpIdx - lfpSzOn;
% keep unshifted time axes, e.g. for comparison:
% tECoG0=tECoG;
% tLFP0=tLFP;

% list of duratons observed in each device:
% dts = change in time [sec]
% dti = change in time [indices]
% differences in these durations tell us how much to adjust
% the time axes
dti_lfp = diff(lfpIdx0);
dts_lfp = dti_lfp/lfpFs;
dti_ecog = diff(ecogIdx);
dts_ecog = dti_ecog/ecogFs;

% at each sync point, jump time ahead with a shift
count_ecog = 0;
count_lfp = 0;
for j = 1:length(dts_lfp)
  disp(['Sync point ' num2str(j)]);
  count_ecog = count_ecog+dti_ecog(j);
  count_lfp = count_lfp+dti_lfp(j);
  shift_j = abs(dts_lfp(j)-dts_ecog(j));
  
  % decide which time axis to shift
  if dts_lfp(j)>dts_ecog(j)
    disp(['ECoG is behind ' num2str(shift_j) ]);
    tECoG(count_ecog:end) = tECoG(count_ecog:end) + shift_j;
  else
    disp(['LFP is behind ' num2str(shift_j)]);
    tLFP(count_lfp:end) = tLFP(count_lfp:end) + shift_j;    
  end
end
toc
% ecogRef0=normalize(ecogRef);
% lfpRef0=normalize(lfpRef(lfpSzOn:lfpSzOff));
% save ~/testSync3.mat -v7.3 patient seizure tLFP tECoG tLFP0 tECoG0 lfpRef0 ecogRef0 lfpSzOn lfpSzOff
% save ~/testSync6.mat -v7.3 patient seizure tLFP tECoG tLFP0 tECoG0 lfpRef0 ecogRef0 lfpSzOn lfpSzOff
% 0;

% Load EEG (if it exists)
if isfield(info, 'EEG')
  tic
  if isfield(info.EEG, 'Labels')
        eegCh = label2index(info.EEG.Labels, ecogProp.Label);
        [dEEG, ~] = openEDF(info.ECoG.EdfFile, eegCh);
        dEEG = dEEG(ecogSzOn : ecogSzOff, :);
        fprintf('EEG loaded\n');
  else
    fprintf('EEG info incomplete\n');
  end
  toc
else
  fprintf('No EEG info found\n');
end

% Now get the original LFP data
tic
lfpProp = NSX_open(info.LFP.Ns5File);
dLFP = openNSx('read', info.LFP.Ns5File, ['c:' num2str(lfpCh(1)) ':' num2str(lfpCh(end))], ...
                ['t:' num2str(lfpSzOn)  ':' num2str(lfpSzOff)], 'precision', 'double');
dLFP = dLFP.Data';                      
fprintf('LFP loaded\n');
toc

% If there's a size mismatch, drop a time point
if length(tLFP)==size(dLFP,2)+1, tLFP(end)=[]; end

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
syncCheck = struct('ecogRef', ecogRef, 'lfpRef', lfpRef, 'ecogIdx', ecogIdx, 'lfpIdx', lfpIdx0, 'lfpSzOn', lfpSzOn, 'lfpSzOff', lfpSzOff);
sz.sync = syncCheck;

% Save data
tic
fprintf('Saving...');
% save([dataPath '/Synced_Multiscale_Mats/' patient '_' seizure '_LFP_ECoG_EEG.mat'], '-v7.3', 'sz');
save(['~/' patient '_' seizure '_LFP_ECoG_EEG.mat'], '-v7.3', 'sz');
fprintf('done.\n');
toc

end
