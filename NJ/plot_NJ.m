%% get betas for all

for i = setdiff(1:32,25)
  i
  fname = ['/media/Shared/GMF/Documents/BostonU/Research/Data/NJ/prot0/10-26-13-prot0-2_pp_c' num2str(i)];
  load(fname,'m');
  save(fname,'m','p');
end
  

%% display P+I+E plots for all cells in a single session

date = '10-26-13'; protocol = 'prot0'; ind = 2;
% date = '10-3-13'; protocol = 'prot1'; ind = 1;
s = Session(protocol,date,ind);
d = s.PPData;
global DO_CONF_INT
DO_CONF_INT = false;

bad_ind = [1 6 10 15 25:27 31:32];
good_ind = setdiff(1:d.N_channels, bad_ind);
% good_ind = [7 9 13 14 16 17 18 20 23 24 29 30];
% good_ind = [9 14 16 17 18 20:24 29 30]; % from 10-3
count = 1;
% %%
% 
% figure('units','normalized','position',[0 0 1 1]);
for response = good_ind
  all_KS(count,:) = m.KS;
% for response = 10
% % for response = 1:d.N_channels
  load([d.Name '_c' num2str(response)]);
  [response, m.KS]
  if m.KS(1)<0.1
    m.plot(d,p);    
    subplot(311), title(num2str(response));
    for n = 1:3, subplot(3,1,n), hold on; end;
    
  else
    
  end
  pause;
  count = count+1;
end

% %% display beta values in the plane + cov. ellipses
% 
% for response = 1:d.N_channels
%   try
%     load([d.Name '_c' num2str(response)]);
%     cov_ellipse(m,p);
%   catch
%     [];
%   end
%   hold on;
% end


%% display P+I+E plots for a variety of cells
global DO_CONF_INT
global PLOT_COLOR

N_plots = 3;
protocol = 'prot1';
DO_CONF_INT = false;

% inhibitory cells
dispCells = { ...
  '5-28-13', 2, 5; '6-13-13', 1, 5; '6-14-13', 2, 6; '7-23-13', 2, 2; ...
  '7-23-13', 2, 2};

% noise-sensitive cells (all)
% dispCells = { ...
%     '6-13-13', 1, 10; '6-14-13', 2, 2; '6-14-13', 2, 8; '7-18-13', 1, 8; ...
%     '7-18-13', 2, 3; '7-23-13', 1, 4; '7-23-13', 2, 3; '7-23-13', 2, 9; ...
%     '7-24-13', 1, 2; '7-24-13', 2, 4; '7-24-13', 2, 5; '10-3-13', 1, 3; ...
%     '10-3-13', 1, 7; '10-3-13', 1, 13; '10-10-13', 1, 4; '10-17-13', 1, 1; ...
%     '10-17-13', 1, 2; '10-17-13', 1, 9; '10-17-13', 1, 10; '10-17-13', 1, 11; ...
%     '10-17-13', 1, 12; '10-17-13', 1, 21; '12-22-13', 1 1; '12-22-13', 1, 6};

% noise-sensitive cells (1)
% dispCells = { ...
%     '6-14-13', 2, 8; '7-18-13', 1, 8; '7-18-13', 2, 3; '7-23-13', 1, 4; ...
%     '7-23-13', 2, 3; '7-24-13', 1, 2; '10-3-13', 1, 7; '10-17-13', 1, 2; ...
%     '10-17-13', 1, 9; '10-17-13', 1, 10; '10-17-13', 1, 11; ...
%     '10-17-13', 1, 12; '10-17-13', 1, 21; '12-22-13', 1, 6};  
% PLOT_COLOR = 'b';

% noise-sensitive cells (2)
% dispCells = { ...
%     '6-13-13', 1, 10; '7-23-13', 2, 9; '7-24-13', 2, 4; '7-24-13', 2, 5; ...
%     '10-3-13', 1, 13; '10-10-13', 1, 4; '12-22-13', 1, 1};
% PLOT_COLOR = 'r';

% noise-sensitive cells (both 1 and 2)
% dispCells = { '6-14-13', 2, 2; '10-3-13', 1, 3};
% PLOT_COLOR = 'm';

% laser-sensitive cells (1)
% dispCells = { ...
%   '6-13-13', 1, 13; '6-18-13', 2, 2; '7-18-13', 1, 1; '7-23-13', 1, 7; ...
%   '7-23-13', 2, 8; '10-3-13', 1, 15; '10-10-13', 1, 5; '10-17-13', 1, 2; ...
%   '10-17-13', 1, 10; '10-17-13', 1, 11; '10-17-13', 1, 13; '12-22-13', 1, 1};
% PLOT_COLOR = 'b';

% laser-sensitive cells (2)
% dispCells = { ...
%   '5-28-13', 1, 1; '7-23-13', 1, 8; '7-24-13', 2, 3; '10-3-13', 1, 5; ...
%   '10-10-13', 1, 11; '10-17-13', 1, 7};
% PLOT_COLOR = 'r';

% laser-sensitive cells (3)
% dispCells = { ...
%   '5-28-13', 1, 11; '6-13-13', 2, 8; '6-14-13', 1, 6; '7-18-13', 2, 5; ...
%   '7-18-13', 2, 7; '7-18-13', 2, 11; '7-24-13', 1, 4; '10-3-13', 1, 6; ...
%   '10-3-13', 1, 10; '10-17-13', 1, 1}; %NOTE: last is also responsive to 1
% PLOT_COLOR = 'g';

% laser-sensitive cells (1 & 2)
% dispCells = { ...
%   '6-11-13', 2, 6; '7-23-13', 1, 6; '10-3-13', 1, 13; '10-10-13', 1, 4; ...
%   '10-17-13', 1, 8; '10-17-13', 1, 9; '10-17-13', 1, 12; '10-17-13', 1, 15; ...
%   '10-17-13', 1, 17; '10-17-13', 1, 19; '12-22-13', 1, 3; '12-22-13', 1, 4};
% PLOT_COLOR = 'm';


% interact cells (all)
% dispCells = {'6-13-13', 2, 6; '6-14-13', 1, 2; '7-18-13', 1, 6; ...
%   '7-18-13', 2, 7; '7-18-13', 2, 11; '7-23-13', 2, 9; '7-23-13', 2, 10; ...
%   '7-24-13', 1, 2; '10-3-13', 1, 7; '10-3-13', 1, 10; '10-3-13', 1, 15; ...
%   '10-3-13', 1, 19; '10-17-13', 1, 9; '10-17-13', 1, 12; '10-17-13', 1, 15};

% interact cells (1)
% dispCells = {'6-13-13', 2, 6; '6-14-13', 1, 2; ...
%   '7-24-13', 1, 2; '10-17-13', 1, 9; '10-17-13', 1, 12; '10-17-13', 1, 15};
% PLOT_COLOR = 'b';

% interact cells (2)
% dispCells = {'7-18-13', 1, 6; '7-18-13', 2, 7; '7-18-13', 2, 11; ...
%   '7-23-13', 2, 9; '7-23-13', 2, 10; '10-3-13', 1, 7; '10-3-13', 1, 10; ...
%   '10-3-13', 1, 15;  '10-3-13', 1, 19};
% PLOT_COLOR = 'r';

load NJ_dt;
% figure('units','normalized','position',[0 0 1 1]);
for i = 1:size(dispCells,1)
  date = dispCells{i,1};
  ind = dispCells{i,2};
  response = dispCells{i,3};
%   s = Session(protocol,date,ind); d = s.PPData;
%   filename = [d.Name '_c' num2str(response)]
  filename = [date '-prot1-' num2str(ind) '_pp_c' num2str(response)];
  load(filename);
  d = pp_data([],(1:length(m.y))*dt);
  p0 = p;
  p0.covariate_names = p0.covariate_names(1:3);
  p0.covariate_channels = p0.covariate_channels(1:3);
  p0.covariate_knots = p0.covariate_knots(1:3);
  p0.covariate_bases = p0.covariate_bases(1:3);
  p0.covariate_ind = p0.covariate_ind(1:3);
  if m.KS(1)<0.1
    m.plot(d,p0);
%     cov_ellipse(m,p);
  else
    fprintf(['Bad model - ' date '_c' num2str(response) '\n']);
  end
  for n = 1:N_plots, subplot(N_plots,1,n), hold on; end;
end
  



%%
