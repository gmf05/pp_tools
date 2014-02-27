pt_name = 'MG49'; sz_name = 'Seizure36';
s = seizure(pt_name,sz_name);
d = s.LFP.PPData;
d = d.downsample(32);
m = pp_model();
p = pp_params();

response = 1;
p.response = response;

T_knots = [0 1];
p = p.add_covar('rate', 0, T_knots, 'indicator');

Q_knots = [1 21:20:101 151:50:501 1001];
p = p.add_covar('self-hist', response, Q_knots, 'spline');

R_knots = [1 21:20:101 151:50:501 1001];
p = p.add_covar('ensemble', [1:response-1, response+1:d.N_channels], R_knots, 'spline');

fit_method = 'glmfit'; noise = [];
% fit_method = 'filt';
% fit_method = 'smooth';
% noise = [0 1e-8 1e-10]; % small, seems to work well
% noise = [0 1e-6 1e-8]; % small, seems to work well
% noise = [0 1e-8 1e-10];
% noise = [0 1e-8 1e-10];
% noise = [0 1e-4 1e-5]

p.fit_method = fit_method;
p.noise = noise;
p.downsample_est = 20;

ms = cell(1,d.N_channels);
for response = 1:d.N_channels
% for response = 1
  response
  p.response = response;
  p.covariate_channels{2} = response;
  p.covariate_channels{3} = [1:response-1, response+1:d.N_channels];
  m = pp_model();
  try 
    m = m.fit(d,p);
    m.X = [];
  catch
    fprintf(['Bad channel\n\n']);
  end 
%   m.plot(d,p); pause;
  ms{response} = m;
end

save([pt_name '_' sz_name '_' p.fit_method],'ms','p');
