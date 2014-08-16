%% load data
patient_name = 'MG49';
seizure_name = 'Seizure45';
data_type = 'LFP';
d1thresh = 0.5; d2thresh = 1;
N = Neuroport(patient_name);

d = get_spikes2(patient_name,seizure_name,data_type,d1thresh,d2thresh);
% d = d.remove_outlier_counts();
dall = d.sub_time(110,160);
dsmall = dall; 
dbig = get_big_spikes(dall,0.1);
dsmall.dn = dsmall.dn - dbig.dn;

%% 
% d = dall;
% d = dbig;
d = dsmall;
d.dn = [d.dn; dbig.dn];

% set knots
dt_ms = round(.001 / d.dt);
T_knots = [0 1]; T_basis = 'indicator';
% T_knots = [0 1]; T_basis = 'spline';
% Q_knots = [] * dt_ms; Q_basis = '';
% R_knots = []; R_basis = '';
Q_knots = [1 30 70 100 200 500] * dt_ms; Q_basis = 'spline';
% Q_knots = [1 300 600 900] * dt_ms; Q_basis = 'spline';
R_knots = [0 5 20]  * dt_ms; R_basis = 'spline'; 
Q = length(Q_knots); R = length(R_knots);

% get list, count of interior electrodes
% use it to count total covariates
chans = cellstr2num(d.labels);
int_elec = N.interior();
N_int = length(int_elec);
N_cov = 1 + 2*(Q+2)*(Q>0) + 4*(R+2)*(R>0); % baseline rate + self-history + spatial
p=[]; ms = cell(1,6);
% cols = {'b','r','k','m','g'};

count=1;
% for response = 1:d.N_channels
for response = int_elec
% for response = 45
  response
%   PLOT_COLOR = cols{count};
  m = pp_model();
  p = pp_params();
  p.response = response;
  c_ind = chans(response);
  
  % first check that it's an interior electrode, then get up/down/etc
  if ismember(c_ind, int_elec)
    [c_up, c_down,c_left,c_right] = N.neighbors(c_ind);
    C_up = find(chans==c_up);
    C_down = find(chans==c_down);
    C_left = find(chans==c_left);
    C_right = find(chans==c_right);
%     D0 = d.sub_data([response,C_up,C_down,C_left,C_right]);
    
    p = p.add_covar('rate',0,T_knots,T_basis); % baseline rate
    if Q>0
      p = p.add_covar('self-history',response,Q_knots,Q_basis);
%       p = p.add_covar('self-history2',response+96,Q_knots,Q_basis);
    end
    
    if R>0
      p = p.add_covar('pop-hist1',C_up,R_knots,R_basis);
      p = p.add_covar('pop-hist2',C_down,R_knots,R_basis);
      p = p.add_covar('pop-hist3',C_left,R_knots,R_basis);
      p = p.add_covar('pop-hist4',C_right,R_knots,R_basis);
    end
    
    m = m.fit(d,p); m, %m.plot(d,p); pause(0.1); %clf;
%     m = m.fit(d,p); m, m.gof(d); pause; clf;
    m.X = []; m.y=[]; ms{count} = m;   
    count=count+1;
  end
end

save([d.name '_chanmodels.mat'],'-v7.3','ms','p');

