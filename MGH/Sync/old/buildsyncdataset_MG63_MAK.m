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

ecogSyncCh = 129; % channel 1 in ECoG carries sync events
lfpSyncCh = 97; % channel 97 in Neuroport LFP carries sync events
[ecogRef, ecogProp] = openEDF(info.ECoG.EdfFile, ecogSyncCh); 
ecogFs = ecogProp.SampleRate(1);
%lfpProp = NSX_open(info.LFP.Ns5File);
%lfpRef = NSX_read(lfpProp, lfpSyncCh, 1, 0, Inf)';
%lfpFs = 1/lfpProp.Period;
d = openNSx('read', info.LFP.Ns5File, ['c:' num2str(lfpSyncCh) ':' num2str(lfpSyncCh)], 'precision', 'double');
lfpRef = d.Data';
lfpFs = d.MetaTags.SamplingFreq;
[ecogIdx, lfpIdx, ecogRealFs] = syncecoglfp(ecogRef, ecogFs, lfpRef, lfpFs);
lfpMaxIdx = length(lfpRef);
%clear ecogRef lfpRef;
fprintf('Synchronized indexes built');

% Load only the ecogCh ECoG channels (there are also EEG data) 
%[dECoG, ecogProp] = openEDF(info.ECoG.EdfFile, ecogCh);
[d, ecogProp] = openEDF(info.ECoG.EdfFile, ecogCh(1));
fclose(ecogProp.FILE.FID);
ecogSzOn = max([1, round((info.StartTime - onset) * ecogFs)]);
ecogSzOff = min([round((info.EndTime + offset) * ecogFs), length(d)]);
% dECoG = dECoG(ecogSzOn : ecogSzOff, :);
ecogRef = ecogRef(ecogSzOn : ecogSzOff);
% fclose(ecogProp.FILE.FID);
dECoG = zeros(ecogSzOff - ecogSzOn + 1, length(ecogCh));
dECoG(:,1) = d(ecogSzOn : ecogSzOff); clear d;
edfFile = info.ECoG.EdfFile;
parfor i = 2:length(ecogCh)
    [di, prop] = openEDF(edfFile, ecogCh(i));
    dECoG(:,i) = di(ecogSzOn : ecogSzOff); 
    di = []; %#ok<NASGU> % force free memory
    fclose(prop.FILE.FID);
end
fprintf(', ECoG loaded');

if isfield(info, 'EEG')
    if isfield(info.EEG, 'Labels')
        eegCh = label2index(info.EEG.Labels, ecogProp.Label);
        [dEEG, ~] = openEDF(info.ECoG.EdfFile, eegCh);
        dEEG = dEEG(ecogSzOn : ecogSzOff, :);
        fprintf(', EEG loaded');
    end
end

% Find the sync events during the sz
idxSz = find(ecogSzOn <= ecogIdx & ecogIdx <= ecogSzOff);

% Use only the time stamps during seizure to set sampling rate.
if length(idxSz) == 1
    if abs(ecogIdx(idxSz)-ecogIdx(idxSz+1)) < abs(ecogIdx(idxSz)-ecogIdx(idxSz-1))
    %Then go to next ecogIdx, it's nearer.
    ecogRealFsSZ = (ecogIdx(idxSz+1) - ecogIdx(idxSz) + 1) / ...
                   ( lfpIdx(idxSz+1) -  lfpIdx(idxSz) + 1) * lfpFs;
    else
    %Else go to previous ecogIdx, it's nearer.
    ecogRealFsSZ = (ecogIdx(idxSz) - ecogIdx(idxSz-1) + 1) / ...
                   ( lfpIdx(idxSz) - lfpIdx(idxSz-1) + 1) * lfpFs;
    end
else
%Else we have more than two indices to consider.
    ecogRealFs0 = zeros(length(idxSz)-1,1);
    for k=1:length(idxSz)-1
        ecogRealFs0(k) = (ecogIdx(idxSz(k+1)) - ecogIdx(idxSz(k)) + 1) / ...
            (lfpIdx(idxSz(k+1)) - lfpIdx(idxSz(k)) + 1) * lfpFs;
    end
    ecogRealFsSZ = mean(ecogRealFs0);
end

% Time difference between ecogSzOn and the first sync event
dtStartSz = (ecogIdx(idxSz(1)) - ecogSzOn) / ecogRealFsSZ;

% Use that to find lfpSzOn
lfpSzOn = round(lfpIdx(idxSz(1)) - dtStartSz * lfpFs);

% Same for the end
dtEndSz = (ecogSzOff - ecogIdx(idxSz(end))) / ecogRealFsSZ;  %Last sync preceeding end of sz.
lfpSzOff = round(lfpIdx(idxSz(end)) + dtEndSz * lfpFs);

if lfpSzOn < 1 || lfpSzOff > lfpMaxIdx
    error('LFP data are truncated compared to ECoG data.');
end

sync_differences = ((ecogIdx(idxSz)-ecogSzOn)/ecogRealFsSZ - (lfpIdx(idxSz)-lfpSzOn)/lfpFs);
indice_mismatch = round(abs(sync_differences)/(1/ecogRealFsSZ));
bad = find(indice_mismatch);
if length(bad) > 0
    fprintf(['Found bad sync points=' num2str(length(bad))])
    %remove_these = ones(length(dECoG),1);
    shift_amount = 0;
    for k=1:length(bad)
        i0 = ecogIdx(idxSz(bad(k))) - ecogSzOn + shift_amount;
        n_to_add = indice_mismatch(bad(k));
        dECoG   = [dECoG([1:i0],:); NaN(n_to_add, size(dECoG,2)); dECoG([i0+1:length(dECoG)],:)];
        ecogRef = [ecogRef([1:i0]); NaN(n_to_add,1); ecogRef([i0+1:length(ecogRef)])];
        
        if isfield(info, 'EEG')
            dEEG = [dEEG([1:i0],:); NaN(n_to_add, size(dEEG,2)); dEEG([i0+1:length(dEEG)],:)];
        end
        
        shift_amount = shift_amount + n_to_add;
        indice_mismatch = indice_mismatch - (n_to_add);
        
        %n_to_remove = indice_mismatch(bad(k))-1;
        %remove_these(i0+1:i0+1+n_to_remove)=0;
        
        %dECoG = dECoG([1:i0-n_to_remove, i0:length(dECoG)], :);
        %ecogRef = ecogRef([1:i0-n_to_remove, i0:length(dECoG)]);
        %indice_mismatch = indice_mismatch - (n_to_remove-1);
    end
end
%keep_these = find(remove_these);
%dECoG = dECoG(keep_these,:);
%ecogRef = ecogRef(keep_these);

% Now get the correct LFP data
%dLFP = NSX_read(lfpProp, lfpCh(1), lfpCh(end), lfpSzOn, lfpSzOff - lfpSzOn + 1, 'p', 'p')';
%fclose(lfpProp.FID);
dLFP = openNSx('read', info.LFP.Ns5File, ['c:' num2str(lfpCh(1)) ':' num2str(lfpCh(end))], ...
                                         ['t:' num2str(lfpSzOn)  ':' num2str(lfpSzOff)], 'precision', 'double');
dLFP = dLFP.Data';
fprintf(', LFP loaded');

ecog = struct('Name', 'ecog', 'File', ecogProp.FileName, 'SamplingRate', ecogRealFs, ...
    'SamplingRateSZ', ecogRealFsSZ, ...
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

syncCheck = struct('ecogRef', ecogRef, 'lfpRef', lfpRef, 'ecogIdx', ecogIdx, 'lfpIdx', lfpIdx);
sz.sync = syncCheck;

fprintf(', saving');
save([dataPath '/Synced_Multiscale_Mats/' patient '_' seizure '_LFP_ECoG_EEG.mat'], '-v7.3', 'sz');
fprintf(', done.\n');

% Plot the syncing results.

% Get LFP ref in the same interval.
%lfpRef = NSX_read(lfpProp, lfpSyncCh, 1, lfpSzOn, lfpSzOff - lfpSzOn + 1, 'p', 'p')';
%fclose(lfpProp.FID);
d = openNSx('read', info.LFP.Ns5File, ['c:' num2str(lfpSyncCh) ':' num2str(lfpSyncCh)], ...
                                         ['t:' num2str(lfpSzOn)  ':' num2str(lfpSzOff)], 'precision', 'double');
lfpRef = d.Data';

ecogFS0 = ecogRealFsSZ;

tECoG = (1:length(ecogRef))/(ecogFS0);
tLFP  = (1:length(lfpRef))/lfpFs;

tsync = tECoG(ecogIdx(idxSz)-ecogSzOn);

for k=1:length(tsync)
    subplot(length(tsync),1,k)
    plot(tECoG, (ecogRef-nanmean(ecogRef))/nanstd(ecogRef), '*')
    hold on
    plot(tLFP, zscore(lfpRef), 'r')
    hold off
    xlim([tsync(k)-0.02, tsync(k)+0.02])
end

set(gcf,'PaperPositionMode','auto');
print('-depsc2', '-tiff', [dataPath '/Synced_Multiscale_Mats/' patient '_' seizure '_syncing.eps']);

end
