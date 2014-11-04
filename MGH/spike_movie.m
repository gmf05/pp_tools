%% get data
% data parameters
% patient_name = 'MG49'; seizure_name = 'Seizure36';
% patient_name = 'MG49'; seizure_name = 'Seizure45';
patient_name = 'MG63'; seizure_name = 'Seizure4';
data_type = 'LFP';
% zthresh = -0.1; dT = 6000;
% zthresh = -0.1; dT = 900;
N = Neuroport(patient_name);

% get raw voltage data
load([DATA '/' patient_name '/' patient_name '_' seizure_name '_LFP_ECoG_EEG.mat']);
dsz = sz.LFP.Data';
tsz = sz.LFP.Time;

% d = get_spikes(patient_name,seizure_name,data_type,0.8,200,0.4,200);

%% plot each spike sweep
close all, figure();
dsFactor = 32;

% % (1) ~78.2: sentinel spike?
% % (2) Omar's movie (~t=750) corresponds to ~t=102.5 in seizure 36
% % sz 36 = [647.5855, 809.6055]

% tmin = 78; tmax = 79;
% tmin = 80; tmax = 81;
tmin = 42.5; tmax = 43.5;
% tmin = 50; tmax = 51;
% tmin = 69; tmax = 70; % -- mg49,s45 beginning of spikes
% tmin = 75; tmax = 76; % 65-75??
% tmin = 102.4; tmax = 104;
% tmin = 125; tmax = 125.5;
% tmin = 121; tmax = 122;
% tmin = 80.2; tmax = 80.8; % mg 63

% voltage data
ti1 = getclosest(tsz,tmin);
ti2 = getclosest(tsz,tmax);
dnorm = -normalize(dsz(:,ti1:ti2));
% dnorm = -zscore(dsz(:,ti1:ti2)')';

d0 = d.sub_time_fast(tmin,tmax).downsample(dsFactor);
% dn0 = cumdownsample(d0.dn,dsFactor);

mov = N.plot(dnorm(:,1:dsFactor:end),[],d0.t,d0.dn);
% mov = N.plot(dnorm(:,1:dsFactor:end),[-4 4],d0.t);
% mytitle = 'test.gif';
mytitle = [patient_name '-' seizure_name '-' num2str(tmin) '-' num2str(tmax) '.gif'];
% movie2gif(mov, mytitle, 'LoopCount', 0, 'DelayTime', 0);



