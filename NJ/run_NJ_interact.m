% indicator functions for stimulus response
% + interaction term
%

% date = '10-10-13'; protocol = 'prot1'; ind = 1;
result = 'interact';
s = Session(protocol,date,ind);

fpref = [date '-' protocol '-' num2str(ind) '_pp_c'];
% files = dir([fpref '*.mat']);
% for c = 1:length(files)
%    load(files(c).name);
%    [c, m.stats.p(end-1:end)']
% end

d = s.PPData;
% d.dn(end+(1:3),:) = 0;
% stim_dur_bins = 0.5/d.dt;
d.dn(d.N_channels+(1:7),:) = 0;
window_bins = floor(0.25/d.dt);

noise_ind = getclosest(d.t,s.noise_TS);
for i = noise_ind    
  d.dn(d.N_channels+1,i-1+(1:window_bins)) = 1;
  d.dn(d.N_channels+2,i-1+window_bins+(1:window_bins)) = 1;
end

laser_ind = getclosest(d.t,s.laser_TS);
for i = laser_ind
  d.dn(d.N_channels+3,i-1+(1:window_bins)) = 1;
  d.dn(d.N_channels+4,i-1+window_bins+(1:window_bins)) = 1;
  d.dn(d.N_channels+5,i-1+2*window_bins+(1:window_bins)) = 1;
end

% interaction terms
d.dn(d.N_channels+6,:) = d.dn(d.N_channels+1,:) .* d.dn(d.N_channels+4,:);
% disp('Paused');
% pause();
% for each noise, if there is laser in the first interval
% set second interval = 1 
ind = noise_ind + window_bins;
for i = ind
  if d.dn(d.N_channels+6,i-10)
    d.dn(d.N_channels+7,i+(1:window_bins)) = 1;
  end
end

p = pp_params();
response = 1;
p.response = response;
T_knots = [0 1];
p = p.add_covar('rate', 0, T_knots, 'indicator');
Q_knots = [1 21:20:101 151:50:501 1001];
p = p.add_covar('self-hist', response, Q_knots, 'spline');
R_knots = [1 21:20:101 151:50:501];
p = p.add_covar('ensemble', [1:response-1, response+1:d.N_channels], R_knots, 'spline');
p = p.add_covar('noise1', d.N_channels+1, 0, 'indicator');
p = p.add_covar('noise2', d.N_channels+2, 0, 'indicator');
p = p.add_covar('laser1', d.N_channels+3, 0, 'indicator');
p = p.add_covar('laser2', d.N_channels+4, 0, 'indicator');
p = p.add_covar('laser3', d.N_channels+5, 0, 'indicator');
p = p.add_covar('interact1',d.N_channels+6, 0, 'indicator');
p = p.add_covar('interact2',d.N_channels+7, 0, 'indicator');

fit_method = 'glmfit'; noise = [];
% fit_method = 'filt';
% fit_method = 'smooth';
% noise = [0 1e-8 1e-10]; % small, seems to work well
% noise = [0 1e-6 1e-8]; % small, seems to work well % noise = [0 1e-8 1e-10]; % noise = [0 1e-8 1e-10]; 
p.fit_method = fit_method;
p.noise = noise;
p.downsample_est = 1000;

% for response = 1
for response = 1:d.N_channels
  p.response = response;
  p.covariate_channels{2} = response;
  p.covariate_channels{3} = [1:response-1 response+1:d.N_channels]; 
  m = pp_model();
%   m = m.make_X(d,p);
  try 
    m = m.fit(d,p);
    m.X = [];
    save([d.Name '_c' num2str(response) '.mat'], 'm','p','result'); clear m;
  catch
    eval(['!echo ' fpref num2str(response) ' >> badchannels']);
  end
end
