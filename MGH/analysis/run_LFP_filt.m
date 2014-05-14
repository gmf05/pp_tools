% pp_tools;
% sz = seizure('MG49','Seizure45');
d = sz.LFP.PPData;
% patient = 'MG49'; seizure = 'Seizure45'; 
% load([DATA_DIR '/' patient '/' patient '_' seizure '_LFP_pp_thresh1']); d = data; clear data
% d = sz.LFP.PPData;
% d = d.downsample(32);
p = pp_params();

response = 44;
p.response = response;

T_knots = [0 1];
p = p.add_covar('rate', 0, T_knots, 'indicator');

Q_knots = [1 21:20:101 151:50:501 1001];
p = p.add_covar('self-hist', response, Q_knots, 'spline');

R_knots = [1 21:20:101 151:50:501];
p = p.add_covar('ensemble', [1:response-1, response+1:d.N_channels], R_knots, 'spline');

% % plane wave dynamics
% R_knots = [0];
% p = p.add_covar('ensemble1', [42], R_knots, 'indicator');
% p = p.add_covar('ensemble2', [43], R_knots, 'indicator');
% p = p.add_covar('ensemble3', [46], R_knots, 'indicator');
% p = p.add_covar('ensemble4', [51], R_knots, 'indicator');



% fit_method = 'glmfit'; noise = [];
fit_method = 'filt';
% fit_method = 'smooth';
noise = [0 1e-6 1e-7]; % gives 10Hz recruitment sig in c1
% noise = [0 1e-6 1e-8]; % 
% noise = [0 1e-6 1e-10]; % previously optimized (needs to be repeated)

p.fit_method = fit_method;
p.noise = noise;
p.downsample_est = 200;

m = pp_model();

bs = cell(1, d.N_channels);
% for response = 1:d.N_channels
for response = 5:10:d.N_channels
  response
  p.response = response; p.covariate_channels{2} = response;
  p.covariate_channels{3} = [1:response-1, response+1:d.N_channels];  
  m = m.fit(d,p);
  m
%   figure, m.plot(d,p); pause; hold off;
%   figure, m.gof(d); pause; close all;
%   bs{response} = m.b;
%   ms{response} = m;

  m.plot(d,p);
  subplot(312); caxis([0,2]);
  subplot(313); caxis([0.8,1.2]);
  pause;
  plot2svg(['newest_filt/mg49_s45_c' num2str(response) '_hinoise.svg']);
  
  subplot(312), ylim([0,200]);
  subplot(313), ylim([0,200]);
  pause;
  plot2svg(['newest_filt/mg49_s45_c' num2str(response) '_hinoise_zoom.svg']);
  close all;
  
end

% % save BU1_S1_static bs ms
% save BU1_S1_static bs

