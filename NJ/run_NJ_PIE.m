% no stimulus response
% just Poisson + Intrinsic + Ensemble effects
%

date = '10-26-13';
protocol = 'prot0';
ind = 2;
s = Session(protocol,date,ind);
d = s.PPData;

fpref = [date '-' protocol '-' num2str(ind) '_pp_c'];
% files = dir([fpref '*.mat']);
% for c = 1:length(files)
%    load(files(c).name);
%    [c, m.stats.p(end-1:end)']
% end

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
p.downsample_est = 1000;

for response = 1:d.N_channels
% for response = 1
  p.response = response;
  p.covariate_channels{2} = response;
  p.covariate_channels{3} = [1:response-1 response+1:d.N_channels]; 
  m = pp_model();
  try 
    m = m.fit(d,p);
    m.X = [];
    save([d.Name '_c' num2str(response) '.mat'], 'm'); clear m;
  catch
    eval(['!echo ' fpref num2str(response) ' >> badchannels']);
  end
end
