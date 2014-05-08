function [ecogIdx, lfpIdx, ecogRealFs] = syncecoglfp_GMF(ecogRef, ecogFs, lfpRef, lfpFs)

% [ECOGIDX, LFPIDX, ECOGREALFS] = SYNCECOGLFP(ECOGREF, ECOGFS, LFPREF,
% LFPFS) uses ECoG and LFP special channels (resp. ECOGREF and LFPREF)
% containing synchronization events to find the corresponding indexes in
% the data (ECOGIDX and LFPIDX), as well as an estimation of the ECoG
% sampling rate ECOGREALFS, assuming that the LFP sampling rate LFPFS is
% correct. ECOGFS is the sampling rate estimated by the recording
% device. Usually ECoG sync events are digitally encoded whereas LFP sync
% events are encoded thanks to pulses. In any cases, a sync event is made
% of a triplet of pulses, the time between pulses within a triplet
% representing a timestamp hh:mm.
% 
% Note: Fs = 500.0056 is perfect for MG49 ecog according to Omar
% This function gives 500.0058 (i.e. 1.4ms difference after 1h) for MG49
% seizure 36
%
% Example
%   [ecogRef, ecogProp] = openEDF(1); % channel 1 in ECoG carries sync events
%   ecogFs = ecogProp.SampleRate(1);
%   lfpProp = NSX_open();
%   lfpRef = NSX_read(lfpProp, 97, 1, 0, Inf)'; % channel 97 in Neuroport LFP
%   lfpFs = 1/lfpProp.Period;
%   [ecogIdx, lfpIdx, ecogRealFs] = syncecoglfp(ecogRef, ecogFs, lfpRef, lfpFs)
% 
% See also: GETSYNCSTAMPS, BUILDSYNCDATASET
%
% Author: Louis-Emmanuel Martinet <louis.emmanuel.martinet@gmail.com>

% Get the timestamps
[ecogStamps, ecogStampsIdx] = getsyncstamps(ecogRef, ecogFs, 0);
z = zeros(length(ecogStamps), 1);
ecogT = datenum([z z z ecogStamps z]);

thresh =  mean(abs(lfpRef)) * 10;
[lfpStamps, lfpStampsIdx] = getsyncstamps(lfpRef, lfpFs, thresh);
z = zeros(length(lfpStamps), 1);
lfpT = datenum([z z z lfpStamps z]);

% Find where is the LFP compared to ECoG
% Be careful that some timestamps are not unique (because every ~30s)

% Which file is starting first between lfp and ecog?
if ecogT(1) < lfpT(1)
    lfpStart = 1;
    ecogStart = find(ecogT == lfpT(lfpStart));
    if length(ecogStart) == 2 && lfpT(lfpStart) ~= lfpT(lfpStart + 1)
        ecogStart = ecogStart(2);
    else 
        ecogStart = ecogStart(1);
    end
else
    ecogStart = 1;
    lfpStart = find(lfpT == ecogT(ecogStart));
    if length(lfpStart) == 2 && ecogT(ecogStart) ~= ecogT(ecogStart + 1)
        lfpStart = lfpStart(2);
    elseif length(lfpStart) == 1 && ecogT(ecogStart) == ecogT(ecogStart + 1)
        ecogStart = 2;
    elseif isempty(lfpStart)
        fprintf('ECoG: %02dh%02d -> %02dh%02d\n', ecogStamps([1 end], :)');
        fprintf('LFP: %02dh%02d -> %02dh%02d\n', lfpStamps([1 end], :)');
        error('No common section in ECoG and LFP data');
    else
        lfpStart = lfpStart(1);
    end
end

% Which file is ending first between lfp and ecog?
if lfpT(end) < ecogT(end)
    lfpEnd = length(lfpT);
    ecogEnd = find(ecogT == lfpT(lfpEnd));
    if length(ecogEnd) == 2 && lfpT(lfpEnd) ~= lfpT(lfpEnd - 1)
        ecogEnd = ecogEnd(1);
    else 
        ecogEnd = ecogEnd(end);
    end
else
    ecogEnd = length(ecogT);
    lfpEnd = find(lfpT == ecogT(ecogEnd));
    if length(lfpEnd) == 2 && ecogT(ecogEnd) ~= ecogT(ecogEnd - 1)
        lfpEnd = lfpEnd(1);
    elseif length(lfpEnd) == 1 && ecogT(ecogEnd) == ecogT(ecogEnd - 1)
        ecogEnd = ecogEnd - 1;
    else
        lfpEnd = lfpEnd(end);
    end
end

% ecogRealFs is the effective sampling rate of all ecog:
ecogIdx = ecogStampsIdx(ecogStart : ecogEnd);
lfpIdx = lfpStampsIdx(lfpStart : lfpEnd);
ecogRealFs = (ecogStampsIdx(ecogEnd) - ecogStampsIdx(ecogStart)) / ...
    (lfpStampsIdx(lfpEnd) - lfpStampsIdx(lfpStart)) * lfpFs;

end

function [stamps, stampsIdx] = getsyncstamps(syncData, fs, thresh)

% [STAMPS, STAMPSIDX] = GETSYNCSTAMPS(SYNCDATA, FS, THRESH) extracts the
% time STAMPS (hour:min and their position STAMPSIDX in the data) of the
% sync event channel SYNCDATA (analog or digital pulse triplets). FS is the
% sampling frequency and THRESH is the threshold to find pulses (e.g.
% mean(abs(SYNCDATA)) * 10 for analog, 0 for digital).
%
% Author: Louis-Emmanuel Martinet <louis.emmanuel.martinet@gmail.com>

% Find the triplet
isAboveThresh = syncData > thresh;
onsetBins = find(isAboveThresh(1:end-1) == 0 & isAboveThresh(2:end) == 1) + 1;
onsetTimes = onsetBins / fs;

% Ensure the first 3 pulses belong to the same triplet
itv = max(diff(onsetTimes(1:4)));
while onsetTimes(3) - onsetTimes(1) > 0.9 * itv
    onsetTimes = onsetTimes(2 : end);
    onsetBins = onsetBins(2 : end);
end
% Same for the end
while onsetTimes(end) - onsetTimes(end - 2) > 0.9 * itv
    onsetTimes = onsetTimes(1 : end - 1);
    onsetBins = onsetBins(1 : end - 1);
end

% Time Conversions
% Hour = ((t2-t1)-100)/50
% Minute = ((t3-t1)-1500)/50
t2t1 = onsetTimes(2 : 3 : end) - onsetTimes(1:3:end);
t3t1 = onsetTimes(3 : 3 : end) - onsetTimes(1:3:end);
hour = round((t2t1 * 1000 - 100) / 50);
min =  round((t3t1 * 1000 - 1500) / 50);
stamps = [hour min];

stampsIdx = onsetBins(diff([onsetTimes(1) - itv ; onsetTimes]) > 0.9 * itv);

end
