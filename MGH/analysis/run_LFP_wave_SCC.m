%--- get data
% % % diary ~/Desktop/MG49_S45_wavemodel_results
% pp_tools;
sz = seizure('MG49','Seizure45');
d = sz.LFP.PPData;
d = sz.LFP.PPData.sub_time(75,115); % seizure 45
% d = sz.LFP.PPData.sub_time(115,150); % seizure 45
% d = sz.LFP.PPData.sub_time(75,110); % seizure 36
% d = sz.LFP.PPData.sub_time(110,145); % seizure 36
% N = Neuroport(sz.Patient);
N = Neuroport('MG49');

%--- set edge/interior electrodes
edge_electrodes = [1:4 6 8 10 14 18 22:2:32 95 65:2:79 81 83 85 88 90 92 93 96];
interior_electrodes = setdiff(1:N.N_electrodes,edge_electrodes);

%--- set some parameters
T_knots = [0 1];
Q_knots = [1 31:10:101 151:50:501 1001];
% R_knots = [0]; % for indicator function
R_knots = [0 50 100 150 250 500]; % for spline function
wave_basis = 'spline';
p = pp_params();
p = p.add_covar('rate', 0, T_knots, 'indicator');
p = p.add_covar('self-hist', -1, Q_knots, 'spline');
p = p.add_covar('ensemble1', -1, R_knots, wave_basis);
p = p.add_covar('ensemble2', -1, R_knots, wave_basis);
p = p.add_covar('ensemble3', -1, R_knots, wave_basis);
p = p.add_covar('ensemble4', -1, R_knots, wave_basis);
bad_electrodes = [14 17 54 63 72 77:79 89]; % added 89 for MG49,Seizure45
% bad_electrodes = [14 54 63 72 77:79 46 56 58 59 68 70 80 87 91];
response_list = [23 25];
% response_list = [11 12 13 23 25];
% response_list = interior_electrodes;
response_list = setdiff(response_list,bad_electrodes);
N_response = length(response_list);
N_covar1 = p.covariate_ind{2}(end);
N_covar2 = length(p.covariate_ind{3}(1):p.covariate_ind{end}(end));

%---initialize arrays
X = zeros(N_response*d.T,N_response*N_covar1+N_covar2);
y = zeros(N_response*d.T,1);

for r = 1:N_response
  response = response_list(r);
  response
  p = pp_params();
  p.response = response;
  p = p.add_covar('rate', 0, T_knots, 'indicator');
  p = p.add_covar('self-hist', response, Q_knots, 'spline');

  %--- plane wave dynamics
  x0 = N.coord(response,1);
  y0 = N.coord(response,2);
  e_up = N.arrayMap(x0,y0+1);
  e_down = N.arrayMap(x0,y0-1);
  e_left = N.arrayMap(x0-1,y0);
  e_right = N.arrayMap(x0+1,y0);

  p = p.add_covar('ensemble1', [e_up], R_knots, wave_basis);
  p = p.add_covar('ensemble2', [e_down], R_knots, wave_basis);
  p = p.add_covar('ensemble3', [e_left], R_knots, wave_basis);
  p = p.add_covar('ensemble4', [e_right], R_knots, wave_basis);
  
  fit_method = 'glmfit'; noise = [];
  % fit_method = 'filt';
  % fit_method = 'smooth';
  % noise = [0 1e-6 1e-7]; % gives 10Hz recruitment sig in c1
  % noise = [0 1e-6 1e-8]; % 
  % noise = [0 1e-6 1e-10]; % previously optimized (needs to be repeated)
  p.fit_method = fit_method;
  p.noise = noise;
  p.downsample_est = 200;

  m = pp_model();
%   m = m.make_X(d,p);
  m = m.fit(d,p);
  m
  X((r-1)*d.T+(p.get_burn_in()+1:d.T),(r-1)*N_covar1+(1:N_covar1)) = m.X(:,1:N_covar1);
%   X((r-1)*d.T+(p.get_burn_in()+1:d.T),end-3:end) = m.X(:,N_covar1+(1:4)); % only for indicator
  X((r-1)*d.T+(p.get_burn_in()+1:d.T),end-N_covar2+1:end) = m.X(:,end-N_covar2+1:end); % only for spline
  X((r-1)*d.T+(1:p.get_burn_in()),(r-1)*N_covar1+1) = 1; % patch rate
  y((r-1)*d.T+(1:d.T),1) = d.dn(response,:)';
  
end

% imagesc(X); 
% plot2svg('designmatrix.svg');
% pause(0.1);

temp = tic();
[b,dev,stats] = glmfit(X,y,'poisson','constant','off');
run_time = toc(temp);
b
run_time

cif = glmval(b,X,'log','constant','off');
LL = sum(log(poisspdf(y,cif)));
% time rescaling & KS test
spike_ind = find(y);           
numISIs = length(spike_ind)-1;
RESCALE_FUNC = 'exp';

if numISIs>2
  z = zeros(1, numISIs);
  for j=1:numISIs                                                
    z(j) = sum(cif(spike_ind(j)+1:spike_ind(j+1)));
  end        
  switch RESCALE_FUNC
    case 'exp'                        
      rs_fn = @(x)(1-exp(-x));
      rs_cdf = @(x)(unifcdf(x,0,1));
    case 'identity'
      rs_fn = @(x)(x);
      rs_cdf = @(x)(expcdf(x,1));
  end
  z = rs_fn(z);

  try
    [eCDF,xCDF] = ecdf(sort(z));
    aCDF = rs_cdf(xCDF);
    ks_stat = max(abs(aCDF-eCDF));
    ks_ci = 1.96/sqrt(numISIs+1);
  catch
    z = [];
    ks_stat = NaN;
    ks_ci = NaN;
  end
else
  z = [];
  ks_stat = NaN;
  ks_ci = NaN;
end
ks_stat

%--- save
% % % file_name = 'MG49_S45_wavemodel_mid_new';
% file_name = 'MG49_S36_wavemodel_mid';
% % % save(file_name,'b','p','cif','LL','z','ks_stat','ks_ci','stats','run_time','response_list');
% % % diary off;

%% plot AR curves

Q_knots = p.covariate_knots{2};
Q = p.covariate_ind{2}(end);
N_response = length(response_list);

figure
for r = 1:N_response
% for r = 58
%   response_list(r)
  ind = (r-1)*Q + (2:Q);
  [t0,y0] = plot_spline(Q_knots,b(ind));
  plot(t0,exp(y0),'b');
  hold on;
%   xlim([0,100]);
%   pause;
end

%% plot population-history curves

% figure

R_knots = p.covariate_knots{3};
N_covar1 = p.covariate_ind{2}(end);
N_response = length(response_list);
start_i = N_covar1*N_response;
N_pop_cov = length(p.covariate_ind{3});

for n = 1:4
  ind = start_i + (1:N_pop_cov);
  subplot(4,1,n); hold on;
  [t0,y0] = plot_spline(R_knots,b(ind));
  plot(t0,exp(y0),'r');
  start_i = start_i + N_pop_cov;
end


%% look at covariance structure

col = stats.covb(:,end-0);
ind = find(col<-.003);
NQ = p.covariate_ind{2}(end);
mod(ind,NQ)

%%
col = stats.covb(:,end-1);
% find largest covariances that AREN'T from first index (i.e. not = 2 mod NQ)
figure, for n = 1:N_response, plot(col((n-1)*NQ+(2:NQ))); hold on; end
ylim([-0.0002,0.0002])

%%

figure
for n = 1:4
  subplot(4,1,n); hold on;
  [t0,y0] = plot_spline(R_knots,b1(p.covariate_ind{2+n}));
  plot(t0,exp(y0),'b');
  [t0,y0] = plot_spline(R_knots,b2(p.covariate_ind{2+n}));
  plot(t0,exp(y0),'r');
end