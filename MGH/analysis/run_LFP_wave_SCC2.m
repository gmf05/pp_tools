%--- get data
diary ~/Desktop/MG49_S45_wavemodel_window_results
% pp_tools;
sz = seizure('MG49','Seizure45');
d0 = sz.LFP.PPData;
N = Neuroport('MG49');

%--- set edge/interior electrodes
edge_electrodes = [1:4 6 8 10 14 18 22:2:32 95 65:2:79 81 83 85 88 90 92 93 96];
interior_electrodes = setdiff(1:N.N_electrodes,edge_electrodes);

%--- set some parameters
T_knots = [0 1];
Q_knots = [1 31:10:101 151:50:501 1001];
% R_knots = [0]; % for indicator function
R_knots = [0 50 100 150 250 500]; % for spline function
wave_basis = 'spline';

p = pp_params();
p = p.add_covar('rate', 0, T_knots, 'indicator');
p = p.add_covar('self-hist', -1, Q_knots, 'spline');
p = p.add_covar('ensemble1', -1, R_knots, wave_basis);
p = p.add_covar('ensemble2', -1, R_knots, wave_basis);
p = p.add_covar('ensemble3', -1, R_knots, wave_basis);
p = p.add_covar('ensemble4', -1, R_knots, wave_basis);
bad_electrodes = [67 89]; % added 89 for MG49,Seizure45
% bad_electrodes = [14 17 54 63 72 77:79 89]; % added 89 for MG49,Seizure45
% bad_electrodes = [14 54 63 72 77:79 46 56 58 59 68 70 80 87 91];
% response_list = [5 7 9 11 12 13 15 16 17 19];
% response_list = [5 7 9 11 13];
% response_list = [9];
response_list = interior_electrodes;
response_list = setdiff(response_list,bad_electrodes);
N_response = length(response_list);
N_covar1 = p.covariate_ind{2}(end);
N_covar2 = length(p.covariate_ind{3}(1):p.covariate_ind{end}(end));

%---initialize arrays
% X = zeros(N_response*d.T,N_response*N_covar1+N_covar2);
% y = zeros(N_response*d.T,1);
window_size = 20; dW = 5;
tmin = 120; tmax = tmin+window_size; % tmax = 140;
tmins = tmin:dW:tmax-window_size;
N_windows = length(tmins);
ms = cell(N_response, N_windows);

file_name = unique_file({sz.Name,'LFP'});

for n = 1:N_windows
  d = d0.sub_time(tmins(n),tmins(n)+window_size);
  for r = 1:N_response
    response = response_list(r);
    response
    p = pp_params();
    p.response = response;
    p = p.add_covar('rate', 0, T_knots, 'indicator');
    p = p.add_covar('self-hist', response, Q_knots, 'spline');

    %--- plane wave dynamics
    x0 = N.coord(response,1);
    y0 = N.coord(response,2);
    e_up = N.arrayMap(x0,y0+1);
    e_down = N.arrayMap(x0,y0-1);
    e_left = N.arrayMap(x0-1,y0);
    e_right = N.arrayMap(x0+1,y0);

    p = p.add_covar('ensemble1', [e_up], R_knots, wave_basis);
    p = p.add_covar('ensemble2', [e_down], R_knots, wave_basis);
    p = p.add_covar('ensemble3', [e_left], R_knots, wave_basis);
    p = p.add_covar('ensemble4', [e_right], R_knots, wave_basis);

    fit_method = 'glmfit'; noise = [];
    % fit_method = 'filt';
    % fit_method = 'smooth';
    % noise = [0 1e-6 1e-7]; % gives 10Hz recruitment sig in c1
    % noise = [0 1e-6 1e-8]; % 
    % noise = [0 1e-6 1e-10]; % previously optimized (needs to be repeated)
    p.fit_method = fit_method;
    p.noise = noise;
    p.downsample_est = 200;

    m = pp_model();
  %   m = m.make_X(d,p);
    try
      m = m.fit(d,p);
      m
      ms{r,n} = m;
    end
    
  %   X((r-1)*d.T+(p.get_burn_in()+1:d.T),(r-1)*N_covar1+(1:N_covar1)) = m.X(:,1:N_covar1);
  %   X((r-1)*d.T+(p.get_burn_in()+1:d.T),end-3:end) = m.X(:,N_covar1+(1:4)); % only for indicator
  %   X((r-1)*d.T+(p.get_burn_in()+1:d.T),end-N_covar2+1:end) = m.X(:,end-N_covar2+1:end); % only for spline
  %   X((r-1)*d.T+(1:p.get_burn_in()),(r-1)*N_covar1+1) = 1; % patch rate
  %   y((r-1)*d.T+(1:d.T),1) d
  end
%   save(file_name,'ms','p','response_list','tmins','tmin','tmax','dW','window_size','N_windows');
end

% % diary off;

%%
D = zeros(N.N_electrodes,4,N_windows);
d = zeros(N.N_electrodes,4);

for n = 1:N_windows
%   figure
%   clf;
  for r = 1:N_response    
    if ~isempty(ms{r,n})      
%       figure(1)
%       for i = 1:4
%         subplot(4,1,i);
%         [t0,y0] = plot_spline(p.covariate_knots{2+i}, ms{r,n}.b(p.covariate_ind{2+i}));
%         d(response_list(r),i) = exp(mean(y0(t0<=450)));
%         if i==1, blah(response_list(r),:) = y0; end
%         y0 = exp(y0);
%         hold on;
%         plot(t0,exp(y0));
%         xlim([0,20]); 
%         ylim([0,40]);
%       end
%       
      figure(2)
      [t0,y0] = plot_spline(p.covariate_knots{2}, ms{r,n}.b(p.covariate_ind{2})); hold on;
%       temp2(r,n) = mean(y0);
      y0 = exp(y0);
      plot(t0,y0); ylim([0,4]); xlim([0,300]);
      title(['MG 49 ' num2str(tmin) '-' num2str(tmax) ' sec']);
%       pause();
      
    end
  end
%   vWobj.writeVideo(getframe(gcf()));
%   pause();
end

%%
% figure, imagesc(blah), caxis([0,2]); colorbar;
% vWobj.close();

% % %%
% % figure, imagesc(temp);
% % 
% % %%
% % for n = 1:N_windows
% %   figure(n); pause;
% % end


% plot direction field using temp
dirv = zeros(N.N_electrodes,2);
dirs = [1 3 2 0];
clear i; % want imaginary i
for r = 1:N_response, z = 0; for j = 1:4, z = z+exp(d(response_list(r),j))*exp(i*dirs(j)*pi/2); end; dirv(response_list(r),:) = [real(z), imag(z)]; end
for j = 1:2, dirv(:,j) = 2*(-1 + 2./(1+exp(-dirv(:,j)))); end

% dirv = [0.5*ones(96,1) pi*2*ones(96,1)];

figure, N.plot_dir(dirv);
% title(['MG 49 ' num2str(tmin) '-' num2str(tmax) ' sec']);
% i=1; dir = 'up';  Ws(response_list) = temp(:,i); N.plot_dir(Ws./max(Ws(Ws>0)),dir);

