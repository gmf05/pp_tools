%% load data
patient_name = 'MG49';
seizure_name = 'Seizure36';
data_type = 'LFP';
thresh = 1;
N = Neuroport(patient_name);

%% get spikes from filtered data
% % 
% % % tmin = 115; tmax = 125;
% % % tmin = 120; tmax = 122;
% % tmin = 120; tmax = 125;
% % % tmin = 100; tmax = 130;
% % % tmin = 80; tmax = 90;
% % min_refract = 0.3*3e4;
% % trange = getclosest(time,tmin):getclosest(time,tmax);
% % time_W = time(trange);
% % Ws = d_post(trange,:)';
% % T = size(Ws,2);
% % labels = str2cell(num2str((1:96)'));
% % 
% % big_dn = zeros(N.N_electrodes,T);
% % small_dn = zeros(N.N_electrodes,T);
% % 
% % for n = 1:N.N_electrodes
% % % for n = [1 10 20]
% %   % get all spikes
% %   [ind,amp] = hilbertspike(Ws(n,:),thresh,1);
% %   % separate into two classes: big and small
% %   [sortAmp,sortI] = sort(amp);
% %   dropI = [];
% %   for j = 2:length(amp)
% %     if min(abs(ind(sortI(1:j-1)) - ind(sortI(j)))) < min_refract
% %       dropI = [dropI sortI(j)];
% %     end
% %   end
% %   small_ind = ind(dropI);
% %   big_ind = ind;
% %   big_ind(dropI) = [];
% % 
% % %   % plot
% % %   plot(time_W,Ws(n,:)); hold on
% % %   plot(time_W(big_ind),Ws(n,big_ind),'rx');
% % %   plot(time_W(small_ind),Ws(n,small_ind),'ks');
% % %   pause(); clf;
% %   
% %   big_dn(n,big_ind) = 1;
% %   small_dn(n,small_ind) = 1;
% % 
% % end
% % 
% % d_big = pp_data(big_dn, time_W, 'name', 'MG49-S36-LFP-big', 'labels', labels);
% % d_small = pp_data(small_dn, time_W, 'name', 'MG49-S36-LFP-small', 'labels', labels);
% % 
% % % d_big = d_big.downsample(4);

load test_wave_model_hi-res

%%

d = d_big;
m = pp_model();
p = pp_params();
N = Neuroport(patient_name);

tcount = 0;
count = 1;

fit_method = 'glmfit'; noise = [];
p.fit_method = fit_method; p.noise = noise;

T_knots = [0 1];
p = p.add_covar('rate', 0, T_knots, 'indicator');

% Q_knots = [1 10:10:100 200:100:500 750 1000];
% Q_knots = [1 20:20:100 200:100:500];
% p = p.add_covar('self-hist', response, Q_knots, basis_fn);

% R_knots = [1 21:20:101 151:50:501];
% p = p.add_covar('ensemble', [1:response-1, response+1:d.N_channels], R_knots, basis_fn);

% R_knots = [0:10:20]; basis_fn = 'spline';
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
%     X(trange,cov_ind) = m.X;
%     y(trange) = d.dn(response,:)';
    
%     m = pp_model();
%     m = m.fit(d,p); m.X=[];
%     m
%     [c_ind,c_up,c_down,c_left,c_right]
%     ps{count} = p;
%     ms{response} = m;

    end
  end
end

% [b,dev,stats] = glmfit(X,y,'poisson','constant','off');
% save wave_model_hi-res b dev stats d p



