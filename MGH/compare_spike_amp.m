% load /projectnb/ecog/Data/MG49/MG49_Seizure45_LFP_filtered
%
figure
thresh = 0.4; spk_style = 'ro';

N_channels = size(d_post,2);
% amps = cell(1,N_channels);

for n = 1:N_channels
  n
% for n = 30
  [spkind,amp] = hilbertspike(d_post(:,n),thresh,1);
  amps2{n} = -d_post(spkind,n);
%   amps{n} = amp;
%   plot(time,d_post(:,n),'b'); hold on;
%   plot(time(spkind),d_post(spkind,n),spk_style);
%   pause; clf;
end
