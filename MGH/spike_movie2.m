% data parameters
patient_name = 'MG49';
seizure_name = 'Seizure45';
data_type = 'LFP';
d1thresh = 0.5;
d2thresh = 1;
% tmin = 120; tmax = 125;
% tmin = 100; tmax = 110;
tmin = 120; tmax = 150;
tdelay = 0.05;

doSave = false;
if doSave
  VW = VideoWriter([patient_name '_' seizure_nam '_new1.avi']);
  VW.open();
end


% get electrode array
N = Neuroport(patient_name);

% get raw voltage data
load([DATA '/' patient_name '/' patient_name '_' seizure_name '_LFP_ECoG_EEG.mat']);
dsz = sz.LFP.Data';
tsz = sz.LFP.Time;
ti1 = getclosest(tsz,tmin);
ti2 = getclosest(tsz,tmax);
tsz = tsz(ti1:ti2);
dsz = dsz(:,ti1:ti2);
mn = min(min(dsz)); mx = max(max(dsz));
dnorm = -normalize(dsz);

% get point process spike data
dpp = get_spikes2(patient_name,seizure_name,data_type,d1thresh,d2thresh);
dpp.labels = str2cell(dpp.labels);
% dpp = dpp.sub_time(tmin,tmax);
% dpp = dpp.remove_outlier_counts();

%% get "big" spikes
dbig = dpp;
dbig.dn = 0*dbig.dn;
count = 1;
for n = 1:dpp.N_channels
  spkind = find(dpp.dn(n,:));
%   mks = normalize(dpp.marks{n}(1,:)); ind = find(mks<0.4); % low amplitude
%   mks = normalize(dpp.marks{n}(3,:)); ind = find(mks>0.5); % high curvature  
%   mks = zscore(dpp.marks{n}(3,:)); ind = find(mks>1); % high curvatur
  mks = dpp.marks{n}(3,:); ind = find(mks>6); % high curvature
  dbig.dn(n,spkind(ind)) = 1;
end

%%


%% plotting

% make axes
close all, figure
ax = cell(1,5);
% left bototm width height
% ax{1} = axes('position',[0.01 0.9 0.98 0.08]);
% ax{2} = axes('position',[0.01 0.8 0.98 0.08]);
% ax{3} = axes('position',[0.01 0.7 0.98 0.08]);
ax{1} = axes('position',[0.01 0.84 0.98 0.12]);
ax{2} = axes('position',[0.01 0.7 0.98 0.12]);
ax{5} = axes('position',[0.03 0.05 0.45 0.58]);
ax{4} = axes('position',[0.51 0.05 0.48 0.58]);


% voltage traces + spikes
% chans = [12 42 80];
chans = [12 42];
Nchans = length(chans);
for n = 1:Nchans
  c = chans(n);
%   figure(1), subplot(Nchans,1,n);
  subplot(ax{n});
  plot(tsz,dsz(c,:)); hold on;
  spkind = find(dpp.dn(c,:)); bigind = find(dbig.dn(c,:));
  plot(tsz(spkind),dsz(c,spkind),'ro','markersize',10,'linewidth',2.5);
  plot(tsz(bigind),dsz(c,bigind),'kx','markersize',10,'linewidth',2.5);
  set(gca,'YTickLabel',[])
  title(['c' num2str(c)]);
  if n<Nchans, set(gca,'XTickLabel',[]);
  else xlabel('time [s]'); end
%   xlim([120,130])  
end
pause;

% get intervals
thresh = 60;
lockout = round(.02/dpp.dt);
I = dbig.spike_trigger(thresh,lockout);
dsort = dpp.sort_mean_time(I);
shiftL = round(0.05/dpp.dt);
shiftR = round(0.15/dpp.dt);
I(:,1) = I(:,1) - shiftL;
I(:,2) = I(:,2) + shiftR;

% for each interval ...
for i = 1:size(I,1)
  ti1 = I(i,1);
  ti2 = I(i,2);
  mn = min(min(dsz(:,ti1:ti2)));
  mx = max(max(dsz(:,ti1:ti2)));
  % draw lines on voltage trace plots to show bounds
  for n = 1:Nchans
    subplot(ax{n});
%     figure(1), subplot(Nchans,1,n)
    h1{n} = plot(tsz([ti1 ti1]), 2*[mn mx], 'r', 'linewidth',2);
    h2{n} = plot(tsz([ti2 ti2]), 2*[mn mx], 'r', 'linewidth',2);
    ylim([1.1*mn 1.1*mx]);
%     xlim([tsz(ti1)-0.05 tsz(ti2)+0.05]);
    xlim([120,125])    
  end
  pause;

  % draw raster plot on re-sorted units
%   figure(2);
  subplot(ax{4});
  cla, dsort.sub_time_fast(ti1:ti2).plot('raster');
  xlim(tsz([ti1 ti2]));
  title('');
  pause;

  % play movie of activity across electrode array
%   figure(2);
  subplot(ax{5});
  title('voltage actiivty on Neuroport array'), %colorbar;
  dsFactor = 30;
  tind = ti1:dsFactor:ti2;
  dn0 = cumdownsample(dpp.dn(:,ti1:ti2),dsFactor);
  N.plot(dnorm(:,tind),[],dpp.t(tind),dn0);
  pause;
    
  %   if doSave
  %     VW.writeVideo(getframe(gcf));
  %   end
  
  % remove lines denoting interval
  for n = 1:Nchans
    delete(h1{n}); delete(h2{n});
  end
end