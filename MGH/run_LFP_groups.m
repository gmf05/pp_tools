pt_name = 'MG49'; sz_name = 'Seizure45';
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

% NOTE: here we are dividing the ensemble into two
% groups based loosely on initial results
% Group 1: left half, less coupled
ens_group1 = [1:5 7 33:45 47:49 51 53 65:88];
% Group 2: right half, more coupled
ens_group2 = [6 8:32 46 50 52 54:64 89:96];

% % % Group 1: left/upper half
% % ens_group1 = [1:2 33 35 37 39:45 47 49 51 53 63 65:93 95:96];
% % % Group 2: right/lower half
% % ens_group2 = [3:32 34 36 38 46 48 50 52 54:62 64 94];

% COVARIATE: ENSEMBLE GROUP 1
R_knots = [1 21:20:101 151:50:501 1001];
p = p.add_covar('ensemble', setdiff(ens_group1,response), R_knots, 'spline');

% COVARIATE: ENSEMBLE GROUP 2
R_knots = [1 21:20:101 151:50:501 1001];
p = p.add_covar('ensemble', setdiff(ens_group2,response), R_knots, 'spline');

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
% for response = 1:d.N_channels
for response = 1
  response
  p.response = response;
  p.covariate_channels{2} = response;
  p.covariate_channels{3} = setdiff(ens_group1,response);
  p.covariate_channels{4} = setdiff(ens_group2,response);  
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

save([pt_name '_' sz_name '_' p.fit_method],'ms','p');
