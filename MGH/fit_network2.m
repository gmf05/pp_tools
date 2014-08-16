patient_name = 'MG49';
seizure_name = 'Seizure45';
data_type = 'LFP';
d1thresh = 0.5;
d2thresh = 1;

N = Neuroport(patient_name);
int_elec = N.interior();
N_int = length(int_elec);

d = get_spikes2(patient_name,seizure_name,data_type,d1thresh,d2thresh);
d.labels = str2cell(d.labels);
% d = d.sub_time(120,150);
% d = d.sub_time(90,120);
d = d.sub_time(60,90);

m = pp_model();
dt_ms = round(.001 / d.dt);

Q_knots = round([1 30:20:100 150:150:750 1000]*dt_ms); Q_basis = 'spline';
R_knots = round([0 5 20]*dt_ms); R_basis = 'spline';
T_knots = [0 1]; T_basis = 'indicator';
Q = length(Q_knots); R = length(R_knots);
N_cov = 1+Q+2+4*(R+2);
chans = cellstr2num(d.labels);
int_elec = N.interior();
ps = [];

% 
X = zeros(N_int*d.T,N_cov);
y = X(:,1);

count = 0;
for response = N.interior()'
% for response = 21:40
% for response = 9
  response
  m = pp_model();
  p = pp_params();
  p.response = response;
%   c_ind = chans(response);
  c_ind = response;
  
  % first check that it's an interior electrode, then get up/down/etc
%   if ismember(c_ind, int_elec)
    [c_up, c_down,c_left,c_right] = N.neighbors(c_ind);
    C_up = find(chans==c_up);
    C_down = find(chans==c_down);
    C_left = find(chans==c_left);
    C_right = find(chans==c_right);
    
    p = p.add_covar('rate',0,T_knots,T_basis); % baseline rate

    if Q>0
      p = p.add_covar('self-history',response,Q_knots,Q_basis);
    end
    
    if R>0
      p = p.add_covar('pop-hist1',C_up,R_knots,R_basis);
      p = p.add_covar('pop-hist2',C_down,R_knots,R_basis);
      p = p.add_covar('pop-hist3',C_left,R_knots,R_basis);
      p = p.add_covar('pop-hist4',C_right,R_knots,R_basis);
    end

    m = m.makeX(d,p);
    trange = count*d.T + (1:d.T);
    X(trange,:) = m.X;
    y(trange) = d.dn(response,:)';
    count = count+1;
    
%     m = m.fit(d,p); m, %m.plot(d,p); pause(); %clf;
%     m = m.fit(d,p); m, m.gof(d); pause; clf;
    
    ps = p; % save parameters
%   end
end
p = ps;

% X0 = X(1:d.T,:); y0 = y(1:d.T);
[b,dev,stats] = glmfit(X,y,'poisson','constant','off');
W = stats.covb;
C = glmval(b,X,'log','constant','off');
m = pp_model(b,W,X,y,C,'log');
m.fit_method = 'glmfit';
save([patient_name '_' seizure_name '_network2'],'-v7.3','m','p');
