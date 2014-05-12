% load /projectnb/ecog/Data/MG49/MG49_Seizure45_LFP_filtered

figure
thresh = 4; spk_style = 'ro';

% N_channels = size(d_post,2);
% for n = 1:N_channels
for n = 30
  spkind = hilbertspike(d_post(:,n),thresh,1);
  plot(time,d_post(:,n),'b'); hold on;
  plot(time(spkind),d_post(spkind,n),spk_style);
%   pause; clf;
end
