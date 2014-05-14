sz = seizure('BU1','Seizure1');
d = sz.ECoG.PPData;
p = pp_params();

response = 1;
p.response = response;

T_knots = [0 1];
p = p.add_covar('rate', 0, T_knots, 'indicator');

Q_knots = [1 21:20:101 151:50:501 1001];
p = p.add_covar('self-hist', response, Q_knots, 'spline');

R_knots = [1 21:20:101 151:50:501];
p = p.add_covar('ensemble', [1:response-1, response+1:d.N_channels], R_knots, 'spline');

fit_method = 'glmfit'; noise = [];
% fit_method = 'filt';
% fit_method = 'smooth';
% noise = [0 1e-8 1e-10]; % small, seems to work well
% noise = [0 1e-6 1e-8]; % small, seems to work well
% noise = [0 1e-8 1e-10];
% noise = [0 1e-8 1e-10]; 

p.fit_method = fit_method;
p.noise = noise;
p.downsample_est = 5;

m = pp_model();

bs = cell(1, d.N_channels);
for response = 1:d.N_channels
  response
  p.response = response; p.covariate_channels{2} = response;
  p.covariate_channels{3} = [1:response-1, response+1:d.N_channels];  
  m = m.fit(d,p);
%   m
%   figure, m.plot(d,p); pause; hold off;
%   figure, m.gof(d); pause; close all;
  bs{response} = m.b;
%   ms{response} = m;
end

% % save BU1_S1_static bs ms
save BU1_S1_static bs

