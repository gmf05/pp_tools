% load both data sets
load([DATA_DIR '/MG49/MG49_Seizure36_LFP_pp_thresh1']);
d36 = obj;
load([DATA_DIR '/MG49/MG49_Seizure45_LFP_pp_thresh1']);
d45 = obj;
clear obj;
%%
% set parameters for a hierarchy of models
response = 1;
N_channels = 96;
T_knots = [0 1];
Q_knots = [1 21:20:101 151:50:501];
R_knots = [0 21:20:101 151:50:501];

p_null = pp_params();
p_ens0 = pp_params();
p_ens1 = pp_params();
p_ens2 = pp_params();

p_null = p_null.add_covar('rate', 0, T_knots, 'indicator');
p_ens0 = p_ens0.add_covar('rate', 0, T_knots, 'indicator');
p_ens1 = p_ens1.add_covar('rate', 0, T_knots, 'indicator');
p_ens2 = p_ens2.add_covar('rate', 0, T_knots, 'indicator');

p_ens0 = p_ens0.add_covar('self-hist', response, Q_knots, 'spline');
p_ens1 = p_ens1.add_covar('self-hist', response, Q_knots, 'spline');
p_ens2 = p_ens2.add_covar('self-hist', response, Q_knots, 'spline');

p_ens1 = p_ens1.add_covar('ensemble', [1:response-1, response+1:N_channels], R_knots, 'spline');
p_ens2 = p_ens2.add_covar('ensemble1', [1:response-1, response+1:N_channels], R_knots, 'spline');

p_ens2 = p_ens2.add_covar('ensemble2', [1:response-1, response+1:N_channels], R_knots, 'spline');

% define "ensemble 1" & "ensemble 2"
ens1 = [1:36 38 42 46:2:52 58:61 63 95];
ens2 = setdiff(1:N_channels,ens1);

response_list = 1:96;
N_response = length(response_list);

% models
m_null = pp_model();
m_ens0 = pp_model();
m_ens1 = pp_model();
m_ens2 = pp_model();

% cells to store estimated models 
ms_null = cell(N_response,2);
ms_ens0 = cell(N_response,2);
ms_ens1 = cell(N_response,2);
ms_ens2 = cell(N_response,2);

% for each electrode, fit model hierarchy
%%
for n = 15:N_response
  n
  try
  response = response_list(n);

  % update parameters
  p_null.response = response;
  p_ens0.response = response;
  p_ens1.response = response;
  p_ens2.response = response;

  p_ens0.covariate_channels{2} = response;
  p_ens1.covariate_channels{2} = response;
  p_ens2.covariate_channels{2} = response;

  p_ens1.covariate_channels{3} = setdiff(1:N_channels,response);

  p_ens2.covariate_channels{3} = setdiff(ens1,response);
  p_ens2.covariate_channels{4} = setdiff(ens2,response);

  m_null = m_null.fit(d36,p_null);
  m_ens0 = m_ens0.fit(d36,p_ens0);
  m_ens1 = m_ens0.fit(d36,p_ens1);
  m_ens2 = m_ens0.fit(d36,p_ens2);
  m_null.X = []; m_ens0.X = []; m_ens1.X = []; m_ens2.X = [];

  ms_null{n,1} = m_null;
  ms_ens0{n,1} = m_ens0;
  ms_ens1{n,1} = m_ens1;
  ms_ens2{n,1} = m_ens2;

  m_null = m_null.fit(d45,p_null);
  m_ens0 = m_ens0.fit(d45,p_ens0);
  m_ens1 = m_ens0.fit(d45,p_ens1);
  m_ens2 = m_ens0.fit(d45,p_ens2);
  m_null.X = []; m_ens0.X = []; m_ens1.X = []; m_ens2.X = [];

  ms_null{n,2} = m_null;
  ms_ens0{n,2} = m_ens0;
  ms_ens1{n,2} = m_ens1;
  ms_ens2{n,2} = m_ens2;
  end
end
save MG49_model_hierarchy ms_null ms_ens0 ms_ens1 ms_ens2 p_null p_ens0 p_ens1 p_ens2

%% plotting

% load MG49_model_hierarchy 
N_channels = 96;
N_response = 96;
DO_CONF_INT = false;
N_plots = 4;

for n = 1:20:N_response
  n
  for i = 1:N_plots, subplot(N_plots,1,i); hold on; end
  
  % plot P+I curves
%   PLOT_COLOR = 'b';
%   ms_ens0{n,1}.plot(d36,p_ens0);
%   PLOT_COLOR = 'r';
%   ms_ens1{n,1}.plot(d36,p_ens0);
%   PLOT_COLOR = 'm';
%   ms_ens2{n,1}.plot(d36,p_ens0);
%   PLOT_COLOR = 'b--';
%   ms_ens0{n,2}.plot(d45,p_ens0);
%   PLOT_COLOR = 'r--';
%   ms_ens1{n,2}.plot(d45,p_ens0);
%   PLOT_COLOR = 'm--';
%   ms_ens2{n,2}.plot(d45,p_ens0);

  % plot ens1 models
%   PLOT_COLOR = 'b';
%   ms_ens1{n,1}.plot(d36,p_ens1);
%   PLOT_COLOR = 'r';
%   ms_ens1{n,2}.plot(d45,p_ens1);
  oldShape = [3 1]; newShape = [1 1]; remapPlots = [0 0 1];

  % plot ens2 models
  PLOT_COLOR = 'b';
  ms_ens2{n,1}.plot(d36,p_ens2);
  PLOT_COLOR = 'r';
  ms_ens2{n,2}.plot(d45,p_ens2);
  oldShape = [4 1]; newShape = [2 1]; remapPlots = [0 0 1 2];
  
%   reshape_subplot(oldShape,newShape,remapPlots);
  pause(); clf;
  
end
%%
oldShape = [3 1]; newShape = [1 1]; remapPlots = [0 0 1]; % ens1
oldShape = [4 1]; newShape = [2 1]; remapPlots = [0 0 1 2]; % ens2
reshape_subplot(oldShape,newShape,remapPlots);

