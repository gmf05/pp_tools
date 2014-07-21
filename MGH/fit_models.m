%% load data
patient_name = 'MG49';
seizure_name = 'Seizure45';
data_type = 'LFP';
d1thresh = 0.5; d2thresh = 1;
N = Neuroport(patient_name);

d = get_spikes2(patient_name,seizure_name,data_type,d1thresh,d2thresh);
d.labels = str2cell(d.labels);
% d = d.remove_outlier_counts();

%% 

% file_name = 'wave_model_s45';
m = pp_model();

T_knots = [0 1]; T_basis = 'indicator';
% T_knots = [0 1]; T_basis = 'spline';
dt_ms = round(.001 / d.dt);
Q_knots = [] * dt_ms; Q_basis = '';
% Q_knots = [1 100:100:500] * dt_ms; Q_basis = 'spline'; 
Q_knots = [1 100:200:900] * dt_ms; Q_basis = 'spline'; 
% Q_knots = [10 200:200:700] * dt_ms; Q_basis = 'spline'; 
R_knots = []; R_basis = '';
% R_knots = [0 5:10:45]  * dt_ms; R_basis = 'spline'; 
% R_knots = [0 10:10:40] * dt_ms; R_basis = 'spline'; 
Q = length(Q_knots); R = length(R_knots);

% get list, count of interior electrodes
chans = cellstr2num(d.labels);
int_elec = N.interior();
N_int = length(int_elec);
N_spatial_cov = 1*(R + (2*isequal(R_basis,'spline')*(R>0)));
N_int_cov = Q+(2*isequal(Q_basis,'spline')*(Q>0));
% N_int_cov = N_int*(Q+(2*isequal(Q_basis,'spline')*(Q>0)));
% Ncov = 1+N_int_cov+N_spatial_cov;
Ncov = N_int + N_spatial_cov;
NT = d.T;

% initialize design matrix, response process
% X = zeros(N_int*NT,Ncov);
% y = zeros(N_int*NT,1);

ps = [];
count = 1;
% for response = 1:d.N_channels
for response = 45
  response
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
%       p = p.add_covar('self-history2',response+1,Q_knots,Q_basis);
    end
    
    if R>0
      p = p.add_covar('pop-hist1',C_up,R_knots,R_basis);
      p = p.add_covar('pop-hist2',C_down,R_knots,R_basis);
      p = p.add_covar('pop-hist3',C_left,R_knots,R_basis);
      p = p.add_covar('pop-hist4',C_right,R_knots,R_basis);
    end

%     m = m.makeX(d,p); % fprintf('Made design matrix\n');
%     trange = (count-1)*NT + (1:NT);
%     cov_ind = [count Ncov+(-N_spatial_cov+1:0)];
%     count = count+1;
%     X(trange,:) = m.X; % all parameters assumed to be = across grid    
%     varind = [count N_int+(1:N_spatial_cov)];
%     X(trange,varind) = m.X;
%     y(trange) = d.dn(response,:)';
    
    m = m.fit(d,p); m, m.plot(d,p); pause; clf;
%     m = m.fit(d,p); m, m.gof(d); pause; clf;
    
    ps = p; % save parameters
  end
end

% %% network model
% [b,dev,stats] = glmfit(X,y,'poisson','constant','off');
% m = pp_model();
% m.X = X; m.y = y; m.b = b; m.W = stats.covb; m.fit_method = 'glmfit';
% m.CIF = glmval(b,X,'log','constant','off');
% m = m.calcGOF();
% figure, m.plot(d,p)
% figure, m.gof(d)
% save(file_name '_m1','-v7.3','m');