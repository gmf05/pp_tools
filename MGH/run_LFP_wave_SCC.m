% pp_tools;
sz = seizure('MG49','Seizure45');
d = sz.LFP.PPData;
% patient = 'MG49'; seizure = 'Seizure45'; 
% load([DATA_DIR '/' patient '/' patient '_' seizure '_LFP_pp_thresh1']); d = data; clear data
% d = sz.LFP.PPData;
% d = d.downsample(32);
p = pp_params();

response = 44;
% response = 46;
p.response = response;

T_knots = [0 1];
p = p.add_covar('rate', 0, T_knots, 'indicator');

Q_knots = [1 21:20:101 151:50:501 1001];
p = p.add_covar('self-hist', response, Q_knots, 'spline');

% R_knots = [1 21:20:101 151:50:501];
% p = p.add_covar('ensemble', [1:response-1, response+1:d.N_channels], R_knots, 'spline');

% plane wave dynamics
R_knots = [0];
p = p.add_covar('ensemble1', [42], R_knots, 'indicator');
p = p.add_covar('ensemble2', [43], R_knots, 'indicator');
p = p.add_covar('ensemble3', [46], R_knots, 'indicator');
p = p.add_covar('ensemble4', [51], R_knots, 'indicator');

% p = p.add_covar('ensemble1', [50], R_knots, 'indicator');
% p = p.add_covar('ensemble2', [44], R_knots, 'indicator');
% p = p.add_covar('ensemble3', [52], R_knots, 'indicator');
% p = p.add_covar('ensemble4', [56], R_knots, 'indicator');



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

bs = cell(1, d.N_channels);
% for response = 1:d.N_channels
for response = response
%   response
%   p.response = response; p.covariate_channels{2} = response;
%   p.covariate_channels{3} = [1:response-1, response+1:d.N_channels];  
  m = m.fit(d,p);
  m
% %   m.plot(d,p);  
end

%% 

N = Neuroport('MG49');
% response_list = [12 27 33 70 31 91 76 49];
response_list = [12 27 33 70];
N_response = length(response_list);
N_covar1 = p.covariate_ind{2}(end);
X = zeros(N_response*d.T,N_response*N_covar1+4);
y = zeros(N_response*d.T,1);

for r = 1:N_response
  
  response = response_list(r);
  p = pp_params();
  p.response = response;
  T_knots = [0 1];
  p = p.add_covar('rate', 0, T_knots, 'indicator');
  Q_knots = [1 21:20:101 151:50:501 1001];
  p = p.add_covar('self-hist', response, Q_knots, 'spline');
  
  % plane wave dynamics
  R_knots = [0];
  x = N.coord(response,1);
  y = N.coord(response,2);
  e_up = N.arrayMap(x,y+1);
  e_down = N.arrayMap(x,y-1);
  e_left = N.arrayMap(x-1,y);
  e_right = N.arrayMap(x+1,y);
  p = p.add_covar('ensemble1', [e_up], R_knots, 'indicator');
  p = p.add_covar('ensemble2', [e_down], R_knots, 'indicator');
  p = p.add_covar('ensemble3', [e_left], R_knots, 'indicator');
  p = p.add_covar('ensemble4', [e_right], R_knots, 'indicator');
  
  fit_method = 'glmfit'; noise = [];
  % fit_method = 'filt';
  % fit_method = 'smooth';
  % noise = [0 1e-6 1e-7]; % gives 10Hz recruitment sig in c1
  % noise = [0 1e-6 1e-8]; % 
  % noise = [0 1e-6 1e-10]; % previously optimized (needs to be repeated)
  p.fit_method = fit_method;
  p.noise = noise;
  p.downsample_est = 200;

  p
  m = m.fit(d,p); m
  X((r-1)*d.T+(p.get_burn_in()+1:d.T),(r-1)*N_covar1+(1:N_covar1)) = m.X(:,1:N_covar1);
  X((r-1)*d.T+(p.get_burn_in()+1:d.T),end-3:end) = m.X(:,N_covar1+(1:4));
    X((r-1)*d.T+(1:p.get_burn_in()),(r-1)*N_covar1+1) = 1; % patch rate
  y((r-1)*d.T+(1:d.T),1) = d.dn(response,:)';
  
end

%%
[b,dev,stats] = glmfit(X,y,'poisson','constant','off');
b
%%
cif = glmval(b,X,'log','constant','off');
LL = sum(log(poisspdf(y,cif)));
% time rescaling & KS test
spike_ind = find(y);           
numISIs = length(spike_ind)-1;

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


%%

%% plot results

% 

% down = 1;
% left = 2;
% right = 3;
% up = 4;
D = 1;
N = Neuroport('MG49');

for i = 1:4
  x = N.coord(p.covariate_channels{i+2},1); y = N.coord(p.covariate_channels{i+2},2);
  % get color
  col = 'b';
  fill([x-D/2 x-D/2 x+D/2 x+D/2],[y-D/2 y+D/2 y+D/2 y-D/2],col); hold on;
end

%%

X0 = zeros(2*d.T,18+22);
X0(1:d.T,1:18) = X(:,1:18);
X0(p.get_burn_in()+1:d.T,1:18) = X(:,1:18);
X0(d.T+p.get_burn_in()+1:end,19:36) = m.X(:,1:18);
X0(p.get_burn_in()+1:d.T,37:40) = X(:,19:22);
X0(d.T+p.get_burn_in()+1:end,37:40) = m.X(:,19:22);
imagesc(X0)

%%
