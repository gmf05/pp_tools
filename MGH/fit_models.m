%% load data
patient_name = 'MG49';
seizure_name = 'Seizure45';
data_type = 'LFP';
thresh = 1;
N = Neuroport(patient_name);

% data_name = [DATA '/' patient_name '/' patient_name '_' seizure_name '_' data_type];
% load([data_name '_filtered']); d_post = d'; time = t;
load wave_data_S45_test

%% get spikes from filtered data

% tmin = 115; tmax = 125;
% tmin = 120; tmax = 122;
% tmin = 120; tmax = 125;
% tmin = 115; tmax = 130;
tmin = 100; tmax = 125;
% tmin = 80; tmax = 90;
min_refract = 0.3*3e4;
trange = getclosest(time,tmin):getclosest(time,tmax);
time_W = time(trange);
Ws = d_post(trange,:)';
T = size(Ws,2);
labels = str2cell(num2str((1:96)'));

big_dn = zeros(N.N_electrodes,T);
small_dn = zeros(N.N_electrodes,T);

for n = 1:N.N_electrodes
% for n = [1 10 20]
  % get all spikes
  [ind,amp] = hilbertspike(Ws(n,:),thresh,1);
  % separate into two classes: big and small
  [sortAmp,sortI] = sort(amp);
  dropI = [];
  for j = 2:length(amp)
    if min(abs(ind(sortI(1:j-1)) - ind(sortI(j)))) < min_refract
      dropI = [dropI sortI(j)];
    end
  end
  small_ind = ind(dropI);
  big_ind = ind;
  big_ind(dropI) = [];

%   % plot
%   plot(time_W,Ws(n,:)); hold on
%   plot(time_W(big_ind),Ws(n,big_ind),'rx');
%   plot(time_W(small_ind),Ws(n,small_ind),'ks');
%   pause(); clf;
  
  big_dn(n,big_ind) = 1;
  small_dn(n,small_ind) = 1;

end

d_big = pp_data(big_dn, time_W, 'name', 'MG49-S45-LFP-big', 'labels', labels);
d_small = pp_data(small_dn, time_W, 'name', 'MG49-S45-LFP-small', 'labels', labels);
d_both = pp_data([small_dn; big_dn], time_W, 'name', 'MG49-S45-LFP-both', 'labels', labels);
d_all = pp_data(small_dn + big_dn, time_W, 'name', 'MG49-S45-LFP', 'labels', labels);
% save -v7.3 wave_data_S45_test  d_big d_small d_both d_all

%% concatenate intervals of interest

d = d_big;

spks = [];
for n = 1:d.N_channels
  spks = [spks d.t(find(d.dn(n,:)))];
end

[counts, times] = hist(spks,d.t(1):0.02:d.t(end));
ind = find(counts>50);
for n = 1:length(ind)
  tOn = times(ind(n)-1);
  tOff = times(ind(n)+1);
  if n==1
    d0 = d.sub_time(tOn,tOff);
  else
    d0 = d0.concat(d.sub_time(tOn,tOff));
  end
end

d0 = d0.reset_time();
d = d0;


%% remove channels with outlying number of spikes

% d = d_small;
d = d_big;
[good_counts,good_channels] = removeoutliers(sum(d.dn,2));
d = d.sub_data(good_channels);

%% small spike models

% 1- what do intrinsic effects give us? good fit, peak ~100ms
% 2- what do add'l spatial effects (small spikes) give us? worse fit, large CI
% 3- what do add'l spatial effects (big spikes) give us?

m = pp_model();
N = Neuroport(patient_name);

T_knots = [0 1]; T_basis = 'indicator';
dt_ms = round(.001 / d.dt);
Q_knots = [] * dt_ms; Q_basis = '';
% Q_knots = [1 100:100:500] * dt_ms; Q_basis = 'spline'; 
% R_knots = []; R_basis = '';
R_knots = [0 5 10] * dt_ms; R_basis = 'spline'; 
Q = length(Q_knots); R = length(R_knots);

% get list, count of interior electrodes
chans = cellstr2num(d.labels);
int_elec = N.interior();
N_int = length(int_elec);
N_spatial_cov = 4*(R + (2*isequal(R_basis,'spline')*(R>0)));
N_int_cov = Q+(2*isequal(Q_basis,'spline')*(Q>0));
% N_int_cov = N_int*(Q+(2*isequal(Q_basis,'spline')*(Q>0)));
Ncov = 1+N_int_cov+N_spatial_cov;
NT = d.T;

% initialize design matrix, response process
X = zeros(N_int*NT,Ncov);
y = zeros(N_int*NT,1);

ps = [];
count = 1;
for response = 1:d.N_channels
% for response = 45
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
    
    try
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

    m = m.makeX(d,p); % fprintf('Made design matrix\n');
    trange = (count-1)*NT + (1:NT);
    cov_ind = [count Ncov+(-N_spatial_cov+1:0)];
    count = count+1;
    X(trange,:) = m.X;
    y(trange) = d.dn(response,:)';
    
% %     m = m.fit(d,p); m, figure, m.plot(d,p); pause; clf;
    
    ps = p; % save parameters
    end
  end
end

%%
[b,dev,stats] = glmfit(X,y,'poisson','constant','off');
m = pp_model();
m.X = X; m.y = y; m.b = b; m.W = stats.covb; m.fit_method = 'glmfit';
m.CIF = glmval(b,X,'log','constant','off');
m = m.calcGOF();
figure, m.plot(d,p)
figure, m.gof(d)



