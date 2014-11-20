%% get data
% data parameters
% patient_name = 'MG49'; seizure_name = 'Seizure36';
% patient_name = 'MG49'; seizure_name = 'Seizure45';
patient_name = 'MG63'; seizure_name = 'Seizure3';
% patient_name = 'MG63'; seizure_name = 'Seizure4';
data_type = 'LFP';
% zthresh = -0.1; dT = 6000;
% zthresh = -0.1; dT = 900;
N = Neuroport(patient_name);

% get raw voltage data
load([DATA '/' patient_name '/' patient_name '_' seizure_name '_LFP_ECoG_EEG.mat']);
dsz = sz.LFP.Data';
tsz = sz.LFP.Time;

% d = get_spikes(patient_name,seizure_name,data_type,0.8,200,0.4,200);

%%
plot(tsz,dsz(1,:)); hold on;
spk = find(d.dn(1,:));
ind = getclosest(tsz,d.t(spk));
plot(tsz(ind),dsz(1,ind),'rx','markersize',8,'linewidth',2);

spk = find(dbig.dn(1,:));
ind = getclosest(tsz,dbig.t(spk));
plot(tsz(ind),dsz(1,ind),'ko','markersize',10,'linewidth',2);

ylim([-14000,4000]);
% x0 = 60; x1 = 80;
x0 = 0; x1 = 20;
for i = 1:5
  xlim([x0,x1])
  pause;  
  print([patient_name '_' num2str(x0) '_' num2str(x1) '.eps'],'-depsc2');
  x0=x0+20; x1=x1+20;
end

%% plot each spike sweep
close all, figure();
dsFactor = 32;

% % (1) ~78.2: sentinel spike?
% % (2) Omar's movie (~t=750) corresponds to ~t=102.5 in seizure 36
% % sz 36 = [647.5855, 809.6055]

% 65.2 - 65.8: early spike mg49
tmin = 65.2; tmax = 65.8;

% voltage data
ti1 = getclosest(tsz,tmin);
ti2 = getclosest(tsz,tmax);
dnorm = -normalize(dsz(:,ti1:ti2));
% dnorm = -zscore(dsz(:,ti1:ti2)')';

d0 = d.sub_time_fast(tmin,tmax).downsample(dsFactor);
% dn0 = cumdownsample(d0.dn,dsFactor);

mov = N.plot(dnorm(:,1:dsFactor:end),[-1 0],d0.t,d0.dn);
% mov = N.plot(dnorm(:,1:dsFactor:end),[-4 4],d0.t);
% mytitle = 'test.gif';
mytitle = [patient_name '-' seizure_name '-' num2str(tmin) '-' num2str(tmax) '.gif'];
% movie2gif(mov, mytitle, 'LoopCount', 0, 'DelayTime', 0);



