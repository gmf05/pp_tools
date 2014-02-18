%% get betas for all 

%% display P+I+E plots for all cells in a single session

% date = '10-26-13'; protocol = 'prot0'; ind = 2;
date = '10-3-13'; protocol = 'prot1'; ind = 1;
s = Session(protocol,date,ind);
d = s.PPData;

% bad_ind = [1 6 10 15 25:27 31:32];
% good_ind = setdiff(1:d.N_channels, bad_ind);
% 
% good_ind = [9 14 16 17 18 20:24 29 30];

% %%
% 
% figure('units','normalized','position',[0 0 1 1]);
% % for response = good_ind
% for response = 10
% % for response = 1:d.N_channels
%   load([d.Name '_c' num2str(response)]);
%   [response, m.KS]
% %   m.plot(d,p);
% %   for n = 1:3, subplot(3,1,n), hold on; end;
% end

% %% display beta values in the plane + cov. ellipses

for response = 1:d.N_channels
  try
    load([d.Name '_c' num2str(response)]);
    cov_ellipse(m,p);
  catch
    [];
  end
  hold on;
end

%% display P+I+E plots for a variety of cells
N_plots = 3;
protocol = 'prot1';

% noise-sensitive cells
% dispCells = { ...
%   '6-11-13', 1, 1; '6-11-13', 1, 4; '6-13-13', 1, 4; '6-13-13', 1, 10; ...
%   '6-14-13', 2, 2; '6-14-13', 2, 8; '7-18-13', 1, 8; '7-18-13', 2, 3; ...
%   '7-23-13', 1, 4; '7-23-13', 2, 3; '7-23-13', 2, 9};
%
% dispCells = { ...
%   '10-3-13', 1, 1; '10-3-13', 1, 3; '10-3-13', 1, 7; '10-10-13', 1, 2; ...
%   '10-10-13', 1, 4; '10-17-13', 1, 1; '10-17-13', 1, 2; '10-10-13', 1, 4; ...
%   '10-17-13', 1, 1; '10-17-13', 1, 2; '10-17-13', 1, 9; '10-17-13', 1, 10; ...
%   '10-17-13', 1, 11; '10-17-13', 1, 12; '10-17-13', 1, 21; '12-22-13', 1, 1; ...
%   '12-22-13', 1, 6};

% noise-sensitive cells
% dispCells = {};


% dispCells = {};
  
figure('units','normalized','position',[0 0 1 1]);
for i = 1:size(dispCells,1)
  date = dispCells{i,1};
  ind = dispCells{i,2};
  response = dispCells{i,3};
  s = Session(protocol,date,ind); d = s.PPData;
  filename = [d.Name '_c' num2str(response)]
  load(filename);
  p0 = p;
  p0.covariate_names = p0.covariate_names(1:3);
  p0.covariate_channels = p0.covariate_channels(1:3);
  p0.covariate_knots = p0.covariate_knots(1:3);
  p0.covariate_bases = p0.covariate_bases(1:3);
  p0.covariate_ind = p0.covariate_ind(1:3);
  m.plot(d,p0);
  for n = 1:N_plots, subplot(N_plots,1,n), hold on; end;
end
  