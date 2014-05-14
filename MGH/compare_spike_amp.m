%% this seems to get BIG spikes

% load ~/Desktop/test.mat % variables: D = data, time = time axis
thresh = 1; spk_style = 'ro'; n = 1;
[spkind,amp] = hilbertspike(D(:,n),thresh,1);

% figure, plot(time,D(:,n)); hold on;
% plot(time(spkind),D(spkind,n),spk_style);

%%-
R = 0.05;
T = 0:0.1:2*pi;
color_RGB = colormap();
mx_amp = max(amp);

d_amp = diff(amp);

for i = 1:length(spkind)-1
  col = 'b';
  if d_amp(i)>0
    fill(time(spkind(i))+R*cos(T), D(spkind(i),n)+R*sin(T), col);
  end
end

% %%
% for i = 1:length(spkind)-1
%   col_ind = ceil(amp(i)/mx_amp*64);
%   col = color_RGB(col_ind,:);  
%   fill(time(spkind(i))+R*cos(T), D(spkind(i),n)+R*sin(T), col);
% end
