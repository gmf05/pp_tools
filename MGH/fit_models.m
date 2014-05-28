%% load data
patient_name = 'MG49';
seizure_name = 'Seizure45';
data_type = 'LFP';
thresh = 1;
N = Neuroport(patient_name);

% data_name = [DATA '/' patient_name '/' patient_name '_' seizure_name '_' data_type];
% load([data_name '_filtered']); d_post = d'; time = t;
% load test_wave_model_hi-res2
% load wave_data_S45_test

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
% save test_wave_model_hi-res2 d_big d_small d_both d_all
% save -v7.3 wave_data_S45_test  d_big d_small d_both d_all

%%

% d = d_big;
% d = d_small; 
m = pp_model();
p = pp_params();
N = Neuroport(patient_name);

tcount = 0;
count = 1;

fit_method = 'glmfit'; noise = [];
p.fit_method = fit_method; p.noise = noise;

T_knots = [0 1];
p = p.add_covar('rate', 0, T_knots, 'indicator');

% Q_knots = [1 100:100:500]; basis_fn = 'spline';
% p = p.add_covar('self-hist', 0, Q_knots, basis_fn);

% R_knots = [0:10:20]; basis_fn = 'spline';
% p = p.add_covar('ensemble', [1:response-1, response+1:d.N_channels], R_knots, basis_fn);

% R_knots = [0 20]; basis_fn = 'spline';
% R_knots = [0 240 480]; basis_fn = 'spline';
R_knots = [0]; basis_fn = 'indicator';

% get list, count of interior electrodes
chans = zeros(1,d.N_channels);
for n = 1:d.N_channels
  chans(n) = str2double(d.labels{n});
end
int_elec = N.interior();
N_int = length(int_elec);
N_spatial_cov = 4*(length(R_knots) + 2*(isequal(basis_fn,'spline')));
Ncov = N_int+N_spatial_cov;
NT = d.T;

% initialize design matrix, response process
X = zeros(N_int*NT,Ncov);
y = zeros(N_int*NT,1);

for response = 1:d.N_channels
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
    
    try      
    p = p.add_covar('rate',0,[0,1],'indicator'); % baseline rate
% %     p = p.add_covar('self-history',response,Q_knots,basis_fn);
    p = p.add_covar('pop-hist1',C_up,R_knots,basis_fn);
    p = p.add_covar('pop-hist2',C_down,R_knots,basis_fn);
    p = p.add_covar('pop-hist3',C_left,R_knots,basis_fn);
    p = p.add_covar('pop-hist4',C_right,R_knots,basis_fn);

% %     m = m.makeX(d,p);
% %     m, pause();
    
% %     NT = size(m.X,1);
% %     trange = (count-1)*NT + (1:NT);
% %     cov_ind = [count Ncov+(-N_spatial_cov+1:0)];
% %     count = count+1;
% %     X(trange,cov_ind) = m.X;
% %     y(trange) = d.dn(response,:)';
    
    m = pp_model();
    m = m.fit(d,p); %m.X=[];
    m, m.b, pause();
% %     ps{count} = p;
% %     ms{response} = m;

    fprintf('Made design matrix\n');
    end
  end
end

[b,dev,stats] = glmfit(X,y,'poisson','constant','off');
% save wave_model_hi-res2 b dev stats d p



