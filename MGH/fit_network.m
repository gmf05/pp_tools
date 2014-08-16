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
d = dsmall;
d.dn = [d.dn; dbig.dn];

% set knots
dt_ms = round(.001 / d.dt);
T_knots = [0 1]; T_basis = 'indicator';
% T_knots = [0 1]; T_basis = 'spline';
% Q_knots = [] * dt_ms; Q_basis = '';
% R_knots = []; R_basis = '';
Q_knots = [1 30 70 100 200 500] * dt_ms; Q_basis = 'spline';
R_knots = [0 5 20]  * dt_ms; R_basis = 'spline'; 
Q = length(Q_knots); R = length(R_knots);

% get list, count of interior electrodes
% use it to count total covariates
chans = cellstr2num(d.labels);
int_elec = N.interior();
N_int = length(int_elec);
N_cov = 1 + 2*(Q+2)*(Q>0) + 4*(R+2)*(R>0); % baseline rate + self-history + spatial

% initialize design matrix, response process
X = zeros(d.T*N_int,N_cov);
y = zeros(d.T*N_int,1);
p=[];
% cols = {'b','r','k','m','g'};

count = 1;
for response = 1:d.N_channels
% for response = int_elec
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
    
    p = p.add_covar('rate',0,T_knots,T_basis); % baseline rate
    if Q>0
      p = p.add_covar('self-history',response,Q_knots,Q_basis);
      p = p.add_covar('self-history2',response+96,Q_knots,Q_basis);
    end
    
    if R>0
      p = p.add_covar('pop-hist1',C_up,R_knots,R_basis);
      p = p.add_covar('pop-hist2',C_down,R_knots,R_basis);
      p = p.add_covar('pop-hist3',C_left,R_knots,R_basis);
      p = p.add_covar('pop-hist4',C_right,R_knots,R_basis);
    end
    
    m = m.makeX(d,p);
%     m = m.fit(d,p); m, %m.plot(d,p); pause(0.1); %clf;
    
    tind = (count-1)*d.T+(1:d.T);
    X(tind,:) = m.X;
    y(tind) = d.dn(response,:)';
    
    count=count+1;
  end
end

[b,dev,stats] = glmfit(X,y,'poisson','constant','off');
W = stats.covb;
C = glmval(b,X,'log','constant','off');
m = pp_model(b,W,X,y,C,'log');
m.fit_method = 'glmfit';
m.X = []

save([d.name '_networkmodel2.mat'],'-v7.3','m','p');
