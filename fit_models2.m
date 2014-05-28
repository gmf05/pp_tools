%% load data
patient_name = 'MG49';
seizure_name = 'Seizure36';
data_type = 'LFP';
thresh = 1;
N = Neuroport(patient_name);
% load test_wave_model_hi-res2

% d = pp_data([d_big.dn; d_small.dn], d_big.t, 'name', [patient_name '-' seizure_name '-' data_type '-all']);
% d = d.downsample(4);

%%

% d = d_big;
% d = d_small; 
m = pp_model();
p = pp_params();
N = Neuroport(patient_name);
fit_method = 'glmfit'; noise = [];
tcount = 0;
count = 1;
T_knots = [0 1];
% Q_knots = [1 100:100:500]; basis_fn = 'spline';
Q_knots = [1 20:20:100 200:100:500]; basis_fn = 'spline';
% R_knots = [0 240 480]; basis_fn = 'spline';
R_knots = [0 10 50]; basis_fn = 'spline';
% R_knots = [0]; basis_fn = 'indicator';

% % % get list, count of interior electrodes
% % chans = zeros(1,d.N_channels);
% % for n = 1:d.N_channels
% %   chans(n) = str2double(d.labels{n});
% % end
% % int_elec = N.interior();
% % N_int = length(int_elec);
% % N_spatial_cov = 4*(length(R_knots) + 2*(isequal(basis_fn,'spline')));
% % Ncov = N_int+N_spatial_cov;
% % NT = d.T;
% % 
% % % initialize design matrix, response process
% % X = zeros(N_int*NT,Ncov);
% % y = zeros(N_int*NT,1);

for response = d.N_channels/2+1:d.N_channels
  response
  m = pp_model();
  p = pp_params();
  p.response = response;
  c_ind = chans(response-d.N_channels/2);
  
  % first check that it's an interior electrode, then get up/down/etc
  if ismember(c_ind, int_elec)    
    [c_up, c_down,c_left,c_right] = N.neighbors(c_ind);
    C_up = find(chans==c_up);
    C_down = find(chans==c_down);
    C_left = find(chans==c_left);
    C_right = find(chans==c_right);
    
    try      
    p = p.add_covar('rate',0,[0,1],'indicator'); % baseline rate
    p = p.add_covar('self-history',response,Q_knots,'spline');
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
    m, pause();
% %     ps{count} = p;
% %     ms{response} = m;

    end
  end
end

% [b,dev,stats] = glmfit(X,y,'poisson','constant','off');
% save wave_model_hi-res2 b dev stats d p



