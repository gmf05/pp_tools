% sz = seizure('MG49','Seizure45');
% d = sz.ECoG.PPData;
load /projectnb/ecog/Data/MG49/MG49_Seizure45_LFP_pp_thresh1
d = obj; clear obj;
% d = d.downsample(32);
p = pp_params();

response = 1;
p.response = response;

T_knots = [0 1];
p = p.add_covar('rate', 0, T_knots, 'indicator');

Q_knots = [1 21:20:101 151:50:501 1001];
p = p.add_covar('self-hist', response, Q_knots, 'spline');

R_knots = [1 21:20:101 151:50:501];
p = p.add_covar('ensemble', [1:response-1, response+1:d.N_channels], R_knots, 'spline');

% fit_method = 'filt';
fit_method = 'smooth';
% noise = [0 1e-4 1e-4]; % larger, ????
% noise = [0 1e-6 1e-8]; % smaller, KS = 0.957, LL=-3.3720e03
noise = [0 1e-10 1e-10]; % optimized KS ~ 0.06-0.08

p.fit_method = fit_method;
p.noise = noise;
p.downsample_est = 50;

m = pp_model();
% bs = cell(1, d.N_channels);
% for response = 1:d.N_channels
for response = 1:10:d.N_channels
% for response = 1
  response
  p.response = response; p.covariate_channels{2} = response;
  p.covariate_channels{3} = [1:response-1, response+1:d.N_channels];  
  try
    m = m.fit(d,p);
    m
    m.X = [];
%   m.plot(d,p); plot2svg(['MG49_Seizure45_smooth_c' num2str(response)]); pause(0.01); clf;
    m.plot(d,p); pause(0.01); saveas(gcf,['MG49_Seizure45_LFP_pp_smooth_c' num2str(response) '.fig']);; 
  end
%   figure, m.gof(d); pause; close all;
  ms{response} = m;
end

save -v7.3 MG49_Seizure45_LFP_pp_smooth ms p
