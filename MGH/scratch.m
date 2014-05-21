%% load data
patient_name = 'MG49';
seizure_name = 'Seizure36';
data_type = 'LFP';
thresh = 1;
N = Neuroport(patient_name);
data_name = [DATA '/' patient_name '/' patient_name '_' seizure_name '_' data_type];
load([data_name '_filtered']); d_post = d'; time = t;

%%

% tmin = 115; tmax = 125;
% tmin = 120; tmax = 122;
tmin = 100; tmax = 130;
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

d_big = pp_data(big_dn, time_W, 'name', 'MG49-S36-LFP-big', 'labels', labels);
d_small = pp_data(small_dn, time_W, 'name', 'MG49-S36-LFP-small', 'labels', labels);

d_big = d_big.downsample(4);

%%

d = d_big;
m = pp_model();
p = pp_params();

T_knots = [0 1];
p = p.add_covar('rate', 0, T_knots, 'indicator');

% Q_knots = [1 10:10:100 200:100:500 750 1000];
% Q_knots = [1 20:20:100 200:100:500];
% p = p.add_covar('self-hist', response, Q_knots, basis_fn);

% R_knots = [1 21:20:101 151:50:501];
% p = p.add_covar('ensemble', [1:response-1, response+1:d.N_channels], R_knots, basis_fn);

fit_method = 'glmfit'; noise = [];
% fit_method = 'filt';
% fit_method = 'smooth';
% noise = [0 1e-6 1e-7]; 
% noise = [0 1e-6 1e-8]; 
% noise = [0 1e-6 1e-10]; % previously optimized (needs to be repeated)

p.fit_method = fit_method;
p.noise = noise;
p.downsample_est = 200;

% R_knots = [0:10:20]; basis_fn = 'spline';
R_knots = [0]; basis_fn = 'indicator';
N = Neuroport(patient_name);
tcount = 0;
count = 1;

N_good_channels = 0;
for response = 1:d.N_channels
  c_ind = str2double(d.labels{response});
  if max(N.coord(c_ind,:))<10 && min(N.coord(c_ind,:))>1
    N_good_channels = N_good_channels + 1;
  end
end
N_spatial_cov = 4*(length(R_knots) + 2*(isequal(basis_fn,'spline')));
Ncov = N_good_channels+N_spatial_cov;
NT = d.T;
X = zeros(N_good_channels*NT,Ncov);
y = zeros(N_good_channels*NT,1);

for response = 1:d.N_channels
  m = pp_model();
  p = pp_params();
  p.response = response;
  c_ind = str2double(d.labels{response});
  
  % first check that it's an interior electrode, then get up/down/etc
  if max(N.coord(c_ind,:))<10 && min(N.coord(c_ind,:))>1
    c_up = N.arrayMap(N.coord(c_ind,1),N.coord(c_ind,2)+1);
    c_down = N.arrayMap(N.coord(c_ind,1),N.coord(c_ind,2)-1);
    c_left = N.arrayMap(N.coord(c_ind,1)-1,N.coord(c_ind,2));
    c_right = N.arrayMap(N.coord(c_ind,1)+1,N.coord(c_ind,2));
    % get indices in labels
    
    clear C_up C_down C_left C_right
    for j = 1:d.N_channels
      if str2double(d.labels{j})==c_up, C_up = j; end
      if str2double(d.labels{j})==c_down, C_down = j; end
      if str2double(d.labels{j})==c_left, C_left = j; end
      if str2double(d.labels{j})==c_right, C_right = j; end
    end
    try
    p = p.add_covar('rate',0,[0,1],'indicator'); % baseline rate
    p = p.add_covar('pop-hist1',C_up,R_knots,basis_fn);
    p = p.add_covar('pop-hist2',C_down,R_knots,basis_fn);
    p = p.add_covar('pop-hist3',C_left,R_knots,basis_fn);
    p = p.add_covar('pop-hist4',C_right,R_knots,basis_fn);

    m = m.makeX(d,p);
    m
    
    NT = size(m.X,1);
    trange = (count-1)*NT + (1:NT);
    cov_ind = [count Ncov+(-N_spatial_cov+1:0)];
    count = count+1;
    X(trange,cov_ind) = m.X;
    y(trange) = d.dn(response,:)';
    
%     m = pp_model();
%     m = m.fit(d,p); m.X=[];
%     m
%     [c_ind,c_up,c_down,c_left,c_right]
%     ps{count} = p;
%     ms{response} = m;

    end
  end
end

[b,dev,stats] = glmfit(X,y,'poisson','constant','off');

save wave_model_big b dev stats d p

% %%
% int_elec = [];
% count = 1;
% for response = 1:d.N_channels
%   c_ind = str2double(d.labels{response});
%   if max(N.coord(c_ind,:))<10 && min(N.coord(c_ind,:))>1
%     int_elec = [int_elec c_ind];
%   end
% end
% %%
% 

