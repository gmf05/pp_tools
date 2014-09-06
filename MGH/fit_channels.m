%% load data
patient_name = 'MG49';
seizure_name = 'Seizure36';
data_type = 'LFP';
d1thresh = 0.5; d2thresh = 1;
N = Neuroport(patient_name);

d = get_spikes2(patient_name,seizure_name,data_type,d1thresh,d2thresh);
% d = d.remove_outlier_counts();
% dall = d.sub_time(115,155);
% dall = d.sub_time(60,110);
dall = d.sub_time(110,135);

%% 
d = dall;
dt_ms = round(.001 / d.dt);
T_knots = [0 1]; T_basis = 'indicator';
Q_knots = [1 30 70 100 200 500] * dt_ms; Q_basis = 'spline';
R_knots = [0 1 5]  * dt_ms; R_basis = 'spline'; 
% R_knots = [0 5 20]  * dt_ms; R_basis = 'spline'; 
Q = length(Q_knots); R = length(R_knots);
int_elec = N.interior();
N_int = length(int_elec);
N_cov = 1 + 1*(Q+2)*(Q>0) + 4*(R+2)*(R>0); % baseline rate + self-history + spatial
p=[];

%% fitting

count=1;
ms = cell(1,N_int);
for response = int_elec
% for response = 45
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
  ps = p;
  try, m = m.fit(d,p); m,
%   m.plot(d,p);
    ms{count} = m;
  end
  count=count+1;
end
p = ps;
save([d.name '_channelmodels3.mat'],'-v7.3','ms','p');
% save([d.name '_channelmodels3_early.mat'],'-v7.3','ms','p');

%%
