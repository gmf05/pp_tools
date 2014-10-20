%% get data
% data parameters
patient_name = 'MG49'; seizure_name = 'Seizure45';
% patient_name = 'MG63'; seizure_name = 'Seizure3';
data_type = 'LFP';
zthresh = -0.1; dT = 6000;
zthresh = -0.1; dT = 6000;
N = Neuroport(patient_name);

% get raw voltage data
load([DATA '/' patient_name '/' patient_name '_' seizure_name '_LFP_ECoG_EEG.mat']);
dsz = sz.LFP.Data';
tsz = sz.LFP.Time;

d = get_spikes(patient_name,seizure_name,data_type,zthresh,dT); % spike data

%% plot each spike sweep
close all, figure();
dsFactor = 32;

% % ~78.2: sentinel spike?
% % tmin = 78;
% % tmax = 79;

% % 750 - 760 ? this is the interval from omar's movie
% % what does this correspond to rel. to seizure onset??
% % answer: t = ~102.5 in seizure 36
% % sz 36 = [647.5855, 809.6055]
% %
% tmin = 102.4;
% tmax = 104;

tmin = 121; tmax = 122;

% tmin = 80.2; tmax = 80.8; % mg 63

% voltage data
ti1 = getclosest(tsz,tmin);
ti2 = getclosest(tsz,tmax);
dnorm = -normalize(dsz(:,ti1:ti2));
% dnorm = -zscore(dsz(:,ti1:ti2)')';

dsFactor2 = 1;
dsFactor3 = dsFactor*dsFactor2;
d0 = d.sub_time(tmin,tmax); % spike data
d0 = d0.downsample(dsFactor);
d0
d0.t(1)
d0.t(end)
1/d0.dt
% dn0 = cumdownsample(d0.dn,dsFactor2);

mov = N.plot(dnorm(:,1:dsFactor:end),[],d0.t,d0.dn);
% mov = N.plot(dnorm(:,1:dsFactor3:end),[-4 4],d0.t);
% mytitle = 'mg49-s45-115-116.gif';
mytitle = [patient_name '-' seizure_name '-' num2str(tmin) '-' num2str(tmax) '.gif'];
movie2gif(mov, mytitle, 'LoopCount', 0, 'DelayTime', 0);



