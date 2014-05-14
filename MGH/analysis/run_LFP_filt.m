% % %% load data
s = seizure('MG49','Seizure45');
d = s.LFP.PPData;
% % d = d.sub_time(20,d.t(end)-20);
% % % d = d.downsample(8);
% % d = d.downsample(32);
% % % d = d.downsample(64);

% set parameters
m = pp_model();
p = pp_params();

response = 1;
p.response = response;

T_knots = [0 1];
p = p.add_covar('rate', 0, T_knots, 'indicator');

Q_knots = [1 21:20:101 151:50:501 1001];
% Q_knots = [1 21:20:101 151:50:501];
p = p.add_covar('self-hist', response, Q_knots, 'spline');

R_knots = [0 21:20:101 151:50:501];
p = p.add_covar('ensemble', [1:response-1, response+1:d.N_channels], R_knots, 'spline');

% fit_method = 'glmfit'; noise = [];
fit_method = 'filt';
% fit_method = 'smooth';
% noise = [0 1e-8 1e-10]; % small, seems to work well
% noise = [0 1e-6 1e-8]; % small, seems to work well
noise = [0 1e-5 1e-7]; % small, seems to work well
% noise = [0 1e-8 1e-10];
% noise = [0 1e-8 1e-10]; 
% noise = [0 1e-5 1e-5];4,

p.fit_method = fit_method;
p.noise = noise;
p.downsample_est = 200;

% ms = cell(1,d.N_channels);
cd(fit_method)
% for response = 1:d.N_channels
for response = 92:d.N_channels
  response
  p.response = response;
  p.covariate_channels{2} = response;
  p.covariate_channels{3} = [1:response-1, response+1:d.N_channels];
  try 
      m = m.fit(d,p);
      m
      m.X = [];
      m.plot(d,p);
      subplot(3,1,2); colorbar; caxis([0,2]);
      subplot(3,1,3); colorbar; caxis([0.8,1.2]);
      plot2svg(['MG49_S45_' fit_method '_' num2str(response), '.svg']); pause(0.1); 
      subplot(3,1,2); ylim([0,150]); subplot(3,1,3); ylim([0,150]);
      plot2svg(['MG49_S45_' fit_method '_zoom_' num2str(response), '.svg']); pause(0.1); clf;
  catch
    m = pp_model();    
    m.KS = [1 1];
  end
  
%   all_KS(response) = m.KS(1);
  
%   ms{response} = m;
end
% cd('../');
% save(['MG49_Seizure45_' fit_method '_KS.mat'], 'all_KS')
