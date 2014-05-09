pat = 'MG49';
sz = 'S36';
doSave = true;
if doSave
  VW = VideoWriter([pat '_' sz '_late.avi']);
  VW.open();
end

%%
% % % % load /media/Shared/GMF/Documents/BostonU/Research/Data/MG49/MG49_Seizure45_LFP_filtered
% % % load /media/Shared/GMF/Documents/BostonU/Research/Data/MG49/MG49_Seizure36_LFP_filtered
% load /projectnb/ecog/Data/MG49/MG49_Seizure36_LFP_filtered
% % % load /media/Shared/GMF/Documents/BostonU/Research/Data/MG49/MG49_Seizure45_LFP_ECoG_EEG
% % % d_post = sz.LFP.Data;
% % % for i = 1:size(d_post,2)
% % %   d_post(:,i) = zscore(d_post(:,i));
% % % end
% % % time = (1:length(d_post))/3e4 - 20;

N = Neuroport('MG49');

tmin = 120; tmax = 122;
% tmin = 115; tmax = 150;
trange = getclosest(time,tmin):getclosest(time,tmax);
time_W = time(trange);
Ws = d_post(trange,:)';
T = size(Ws,2);

spike_dn = zeros(N.N_electrodes,T);
for n = 1:N.N_electrodes
  ind = hilbertspike(Ws(n,:),1,1);
  spike_dn(n,ind) = 1;
end
% 
% Ws = -Ws;
% 
% % %%
% % T = size(d_post,1);
% % spike_dn = zeros(N.N_electrodes,T);
% % for n = 1:N.N_electrodes
% % % for n = 1
% %   n
% %   ind = hilbertspike(d_post(:,n),1,1);
% %   spike_dn(n,ind) = 1;
% % end
% % %%

%%
for n = 1:N.N_electrodes
% for n = [41 79 82]
  plot(time,d_post(:,n)); hold on;
  spike_ind = find(spike_dn(n,:));
  plot(time_W(spike_ind), Ws(n,spike_ind),'ro');
  pause; hold off;
end


%%
% mn = min(min(Ws));
% mx = max(max(Ws));
% cax = [mn mx];
% cax = [-4 4];
cax = [6 -6];
C = cax(2)-cax(1);
R = 0.5;
theta = 0:0.01:2*pi;

%  figure('units','normalized','position',[0 0 1 1]);
clf, set(gcf,'units','normalized','position',[0 0 1 1]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
colormap('default'), color_RGB = colormap();


% tmn = 84; tmx = 86.5;
tmn = time_W(1); tmx = time_W(end);
response_list = [33, 42, 84];
% response_list = [41 76 82];
for i = 1:3
  r = response_list(i);
  spike_ind = find(spike_dn(r,:));
  subplot(3,2,2*i); plot(time_W,Ws(r,:)); hold on, plot(time_W(spike_ind), Ws(r,spike_ind),'rx');  xlim([tmn,tmx]); ylim([-8,8]);
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
  pause(0.002);
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
