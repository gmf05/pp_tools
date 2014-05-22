patient_name = 'MG49';
seizure_name = 'Seizure36';
data_type = 'LFP';
thresh = 1;
N = Neuroport(patient_name);
data_name = [DATA '/' patient_name '/' patient_name '_' seizure_name '_' data_type];
load([data_name '_filtered']); d_post = d'; time = t;

doSave = false;
if doSave
  VW = VideoWriter([patient_name '_' seizure_nam '_late.avi']);
  VW.open();
end


%%
% tmin = 115; tmax = 125;
% tmin = 120; tmax = 122;
tmin = 120; tmax = 125;
% tmin = 80; tmax = 90;
min_refract = 0.3*3e4;
trange = getclosest(time,tmin):getclosest(time,tmax);
time_W = time(trange);
Ws = d_post(trange,:)';
T = size(Ws,2);
amps = cell(1,N.N_electrodes);
spike_dn = zeros(N.N_electrodes,T);
for n = 1:N.N_electrodes
% for n = [1 10 20]
  [ind,amp] = hilbertspike(Ws(n,:),thresh,1);
%   i = find(diff(amp)<0)+1;
%   i = find(amp<mean(amp));

  [sortAmp,sortI] = sort(amp);
  dropI = [];
  for j = 2:length(amp)
    if min(abs(ind(sortI(1:j-1)) - ind(sortI(j)))) < min_refract
      dropI = [dropI sortI(j)];
    end
  end
  ind(dropI) = [];

  spike_dn(n,ind) = 1;
%   amps{n} = amp;
end

% d = pp_data(spike_dn,time_W);

%%
% i = find(diff(amp)<0);
% i = find(amp<mean(amp));

figure
plot(time_W,Ws(n,:)); hold on
plot(time_W(ind),Ws(n,ind),'rx');
% plot(time_W(ind(i)),Ws(n,ind(i)),'ks');

% spike_dn = zeros(N.N_electrodes,T);
% spike_dn(1,ind(i)) = 1;

%%
% T = size(d_post,1);
% spike_dn = zeros(N.N_electrodes,T);
% for n = 1:N.N_electrodes
% % for n = 1
%   n
% %   ind = hilbertspike(d_post(:,n),1,1);
%   spike_dn(n,ind) = 1;
% end
% %%

% %%
% for n = 1:N.N_electrodes
% % for n = [41 79 82]
%   plot(time,d_post(:,n)); hold on;
%   spike_ind = find(spike_dn(n,:));
%   plot(time_W(spike_ind), Ws(n,spike_ind),'ro');
%   pause; hold off;
% end

% %%
% mn = min(min(Ws));
% mx = max(max(Ws));
% cax = [mn mx];
% cax = [-4 4];
cax = [4 -6];
C = cax(2)-cax(1);
R = 0.5;
theta = 0:0.01:2*pi;

%  figure('units','normalized','position',[0 0 1 1]);
clf, set(gcf,'units','normalized','position',[0 0 1 1]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
colormap('default'), color_RGB = colormap();


% tmn = 105; tmx = 110;
tmn = time_W(1); tmx = time_W(end);
% response_list = [33, 42, 84];
response_list = [1 10 20];
% response_list = [41 76 82];
for i = 1:3
  r = response_list(i);
  subplot(3,2,2*i); plot(time_W,Ws(r,:));
  spike_ind = find(spike_dn(r,:));
  hold on, plot(time_W(spike_ind), Ws(r,spike_ind),'rx');  xlim([tmn,tmx]); ylim([-8,8]);
end

dW = 1;
dT = 30;
for t = dT+1:dT:T
  subplot(1,2,1)
  for n = 1:N.N_electrodes
    x = N.coord(n,1);
    y = N.coord(n,2);
    col_ind = round((Ws(n,t)-cax(1))/C*63)+1;
    col_ind = min(col_ind,64); col_ind = max(col_ind,1);
    col = color_RGB(col_ind,:);
    fill([x-R x-R x+R x+R],[y-R y+R y+R y-R], col); hold on;
%     text(x,y,num2str(n),'fontsize',22);
  end
  
  for r = response_list, text(N.coord(r,1),N.coord(r,2),num2str(r),'fontsize',20); end;
  
  spike_cells = find(sum(spike_dn(:,t-dT:t)'));
  if ~isempty(spike_cells)
    for n = spike_cells, fill(N.coord(n,1)+R/2*cos(theta),N.coord(n,2)+R/2*sin(theta),'w'); end
  end
  
  axis xy;
  
  subplot(3,2,2), hold on; time_bar = plot(time_W(t)*ones(1,2),[-8,8],'r','LineWidth',2);

  subplot(1,2,1);
  title(num2str(time_W(t)));
  pause(0.005);
  hold off;
  
  if doSave
    VW.writeVideo(getframe(gcf));
  end
  delete(time_bar)
end

if doSave
  close(gcf)
  VW.close();
end
