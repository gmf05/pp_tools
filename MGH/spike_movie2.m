% data parameters
patient_name = 'MG49';
seizure_name = 'Seizure45';
data_type = 'LFP';
tmin = 120; tmax = 130;

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
dnorm = normalize(dsz);

% get point process spike data
dpp = get_spikes2(patient_name,seizure_name,data_type,0.5);
dpp = dpp.sub_time(tmin,tmax);

%% get "big" spikes
dn = 0*dpp.dn;
count = 1;
for n = 1:dpp.N_channels
  spkind = find(dpp.dn(n,:));
  mks = normalize(dpp.marks{n}(1,:)); ind = find(mks<0.5);
  dn(n,spkind(ind)) = 1;
end
dbig = dpp;
dbig.dn = dn;

%% plotting

% figure

% voltage traces + spikes
chans = [1 2 3];
Nchans = length(chans)
for n = 1:Nchans
  c = chans(n);
  figure(1), subplot(Nchans,1,n);
  plot(tsz,dsz(c,:)); hold on;
  spkind = find(dpp.dn(c,:)); bigind = find(dbig.dn(c,:));
  plot(tsz(spkind),dsz(c,spkind),'ro','markersize',10,'linewidth',2.5);
  plot(tsz(bigind),dsz(c,bigind),'kx','markersize',10,'linewidth',2.5);
  set(gca,'YTickLabel',[])
  title(['c' num2str(c)]);
end

% get intervals
thresh = 65;
lockout = round(.02/d.dt);
I = dbig.spike_trigger(thresh,lockout);
dsort = dpp.sort_mean_time(I);
shiftL = round(0.2/dpp.dt);
shiftR = round(0.3/dpp.dt);
I(:,1) = I(:,1) - shiftL;
I(:,2) = I(:,2) + shiftR;

% for each interval ...
for i = 1:size(I,1)
  ti1 = I(i,1);
  ti2 = I(i,2);
  
  % draw lines on voltage trace plots to show bounds
  for n = 1:Nchans
    figure(1), subplot(Nchans,1,n);
    c = chans(n);
    h1{n} = plot(tsz([ti1 ti1]), [mn mx], 'r', 'linewidth',2);
    h2{n} = plot(tsz([ti2 ti2]), [mn mx], 'r', 'linewidth',2);
  end
  pause;

  % draw raster plot on re-sorted units
  % get axis
  figure(2), clf, dsort   .sub_time_fast(ti1:ti2).reset_time().plot('raster');
  pause;

  % play spike movie (variable pause delay)
  figure(3)    
  for j = ti1:10:ti2, N.plot(dnorm(:,j)); pause(tdelay); end
  pause;
    
  %   if doSave
  %     VW.writeVideo(getframe(gcf));
  %   end
    
  % remove lines denoting interval
  for n = 1:Nchans
    delete(h1{n}); delete(h2{n});
  end
end