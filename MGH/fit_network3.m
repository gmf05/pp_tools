patient_name = 'MG49';
seizure_name = 'Seizure45';
data_type = 'LFP';
d1thresh = 0.5;
d2thresh = 1;

N = Neuroport(patient_name);
int_elec = N.interior();
N_int = length(int_elec);

d = get_spikes2(patient_name,seizure_name,data_type,d1thresh,d2thresh);
d = d.sub_time(120,150);
d = get_big_spikes(d,0.1);

m = pp_model();
dt_ms = round(.001 / d.dt);

% Q_knots = round([1 30:20:100 150:150:750 1000]*dt_ms); Q_basis = 'spline';
Q_knots = round([1 200:200:1000]*dt_ms); Q_basis = 'spline';
R_knots = round([0 5 20]*dt_ms); R_basis = 'spline';
% Q_knots = []; Q_basis = 'spline';
% R_knots = []; R_basis = 'spline';
T_knots = [0 1]; T_basis = 'indicator';
Q = length(Q_knots); R = length(R_knots);
N_cov = 1+(Q+2)*(Q>0)+4*(R+2)*(R>0);
int_elec = N.interior();
ps = [];

% 
X = zeros(N_int*d.T,N_cov);
y = X(:,1);

count = 0;
for response = N.interior()'
% for response = 9
  response
  m = pp_model();
  p = pp_params();
  p.response = response;
  [c_up, c_down,c_left,c_right] = N.neighbors(response);
  
  p = p.add_covar('rate',0,T_knots,T_basis); % baseline rate
  if Q>0
    p = p.add_covar('self-history',response,Q_knots,Q_basis);
  end

  if R>0
    p = p.add_covar('pop-hist1',c_up,R_knots,R_basis);
    p = p.add_covar('pop-hist2',c_down,R_knots,R_basis);
    p = p.add_covar('pop-hist3',c_left,R_knots,R_basis);
    p = p.add_covar('pop-hist4',c_right,R_knots,R_basis);
  end

  m = m.makeX(d,p);
  trange = count*d.T + (1:d.T);
  X(trange,:) = m.X;
  y(trange) = d.dn(response,:)';
  count = count+1;
    
%     m = m.fit(d,p); m, %m.plot(d,p); pause(); %clf;
%     m = m.fit(d,p); m, m.gof(d); pause; clf;
    
    ps = p; % save parameters
end
p = ps;

% X0 = X(1:d.T,:); y0 = y(1:d.T);
[b,dev,stats] = glmfit(X,y,'poisson','constant','off');
W = stats.covb;
C = glmval(b,X,'log','constant','off');
m = pp_model(b,W,X,y,C,'log');
m.fit_method = 'glmfit';
m.X = []
save([patient_name '_' seizure_name '_bigspikes_network2'],'-v7.3','m','p');
