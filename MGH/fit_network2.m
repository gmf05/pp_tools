%% load data
patient_name = 'MG49';
% seizure_name = 'Seizure45';
seizure_name = 'Seizure36';
data_type = 'LFP';
d1thresh = 0.5; d2thresh = 1;
N = Neuroport(patient_name);

d = get_spikes2(patient_name,seizure_name,data_type,d1thresh,d2thresh);
% d = d.remove_outlier_counts();
% dall = d.sub_time(115,155);
% dall = d.sub_time(60,115);
% dall = d.sub_time(110,135);
dall = d.sub_time(60,110);

%% 
d = dall;
dt_ms = round(.001 / d.dt);
T_knots = [0 1]; T_basis = 'indicator';
Q_knots = [1 30 70 100 200 500] * dt_ms; Q_basis = 'spline';
% R_knots = [0 5 20]  * dt_ms; R_basis = 'spline'; 
R_knots = [0 1 5]  * dt_ms; R_basis = 'spline'; 
Q = length(Q_knots); R = length(R_knots);
int_elec = N.interior();
N_int = length(int_elec);
N_cov = 1 + 1*(Q+2)*(Q>0) + 4*(R+2)*(R>0); % baseline rate + self-history + spatial
X = zeros(d.T*N_int,N_cov);
y = zeros(d.T*N_int,1);
p=[];

%% fitting
count = 0;
% for response = 1:d.N_channels
for response = int_elec
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
  tind = count*d.T+(1:d.T);
  X(tind,:) = m.X;
  y(tind) = d.dn(response,:)';

  count=count+1;
end

[b,dev,stats] = glmfit(X,y,'poisson','constant','off');
W = stats.covb;
C = glmval(b,X,'log','constant','off');
m = pp_model(b,W,X,y,C,'log');
m.fit_method = 'glmfit';
m.X = []

% save([d.name '_networkmodel2_allspikes.mat'],'-v7.3','m','p');
save([d.name '_networkmodel2_allspikes_early.mat'],'-v7.3','m','p');
