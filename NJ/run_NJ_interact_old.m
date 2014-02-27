% indicator functions for stimulus response
% + interaction term
%

% date = '10-17-13'; protocol = 'prot1'; ind = 1;
% result = 'interact';
% s = Session(protocol,date,ind);

fpref = [date '-' protocol '-' num2str(ind) '_pp_c'];
% files = dir([fpref '*.mat']);
% for c = 1:length(files)
%    load(files(c).name);
%    [c, m.stats.p(end-1:end)']
% end

d = s.PPData;
d.dn(end+(1:3),:) = 0;
stim_dur_bins = 0.5/d.dt;

noise_ind = getclosest(d.t,s.noise_TS);
for i = noise_ind    
  d.dn(end-2,i+(0:stim_dur_bins)) = 1;
end

laser_ind = getclosest(d.t,s.laser_TS);
for i = laser_ind
  d.dn(end-1,i+(0:stim_dur_bins)) = 1;
end

% interaction term
d.dn(end,:) = d.dn(end-1,:).*d.dn(end-2,:);

p = pp_params();
response = 1;
p.response = response;
T_knots = [0 1];
p = p.add_covar('rate', 0, T_knots, 'indicator');
Q_knots = [1 21:20:101 151:50:501 1001];
p = p.add_covar('self-hist', response, Q_knots, 'spline');
R_knots = [1 21:20:101 151:50:501];
p = p.add_covar('ensemble', [1:response-1, response+1:d.N_channels], R_knots, 'spline');
p = p.add_covar('noise', d.N_channels+1, 0, 'indicator');
p = p.add_covar('laser', d.N_channels+2, 0, 'indicator');
p = p.add_covar('interact',d.N_channels+3, 0, 'indicator');

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
  p.response = response;
  p.covariate_channels{2} = response;
  p.covariate_channels{3} = [1:response-1 response+1:d.N_channels]; 
  m = pp_model();
  try 
    m = m.fit(d,p);
    m.X = [];
    save([d.Name '_c' num2str(response) '.mat'], 'm','p','result'); clear m;
  catch
    eval(['!echo ' fpref num2str(response) ' >> badchannels']);
  end
end
