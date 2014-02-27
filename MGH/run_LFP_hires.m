s = seizure('MG49','Seizure45');
d = s.LFP.PPData;
m = pp_model();
p = pp_params();

response = 1;
p.response = response;

T_knots = [0 1];
p = p.add_covar('rate', 0, T_knots, 'indicator');

Q_knots = [1 101:200:1001 1501:500:5001 10001];
p = p.add_covar('self-hist', response, Q_knots, 'spline');

R_knots = [1 101:200:1001 1501:500:5001 10001];
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

ms = cell(1,d.N_channels);
for response = 1:d.N_channels
  response
  p.response = response;
  p.covariate_channels{3} = [1:response-1, response+1:d.N_channels];
  m = pp_model();
  try 
    m = m.fit(d,p);
    m.X = [];
    m
    m.b
  catch
    fprintf(['Bad channel\n\n']);
  end 
%   m.plot(d,p); pause;
  ms{response} = m;
end

save MG49_Seizure45_LFP_static_hires ms p
