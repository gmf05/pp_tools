% load MG49_model_hierarchy
sz = seizure('MG49','Seizure45');
d = sz.LFP.PPData;
sz_ind = 2;

%% poisson model
sim_dn = zeros(d.N_channels,d.T);
for i = 1:d.N_channels
  if ~isempty(ms_null{i,sz_ind})
    sim_dn(i,:) = boolean(poissrnd(exp(ms_null{i,sz_ind}.b),[1,d.T]));  
  end  
end
sd = pp_data(sim_dn,d.t);

%% more simulated data: delay network

% set baseline rate
T_knots = [0 1];
T_b = log(30*d.dt);

% express delays as [knot,b] combos
R_knots = [1 5 10 15];
R_b = [-10 -3 0 1 0 0]';

% declare parameters
p_est = pp_params();
p_est.response = 1;
p_est = p_est.add_covar('rate',0,T_knots,'indicator'); % flat baseline rate
p_est = p_est.add_covar('ens',2,R_knots,'spline');

% plot pair (later, sequence) of delay effect curves
[t,y] = plot_spline(R_knots,R_b);
plot(t,y);

% simulate data
% sim_dn = zeros(d.N_channels,d.T);
sim_dn = zeros(2,d.T);

% simulate poisson channel
sim_dn(2,:) = boolean(poissrnd(20*d.dt,[1,d.T]));
sd = pp_data(sim_dn,d.t);
burn_in = p_est.get_burn_in();

m_true = pp_model();
m_true = m_true.make_X(sd,p_est);
m_true.b = [T_b; R_b];
sim_dn(1,:) = poissrnd(exp(m_true.X*m_true.b));
sd = pp_data(sim_dn,d.t);

% estimate delay effect curves
m_est = pp_model();
m_est = m_est.fit(sd,p_est);
R_b_est = m_est.b(p_est.covariate_ind{2});

% compare estimate to ground truth
[t,yhat] = plot_spline(R_knots,R_b_est);
hold on; plot(t,yhat,'r');

%% more simulated data: delay network 2 nodes

N_nodes = 3;

% set baseline rate
T_knots = [0 1];
T_b = log(30*d.dt);

% express delays as [knot,b] combos
R_knots = [1 5 10 15];
R_b = cell(1,N_nodes-1);
R_b{1} = [-10 -3 0 1 0 0]';
R_b{2} = [-10 -3 0 0 1 0]';

% declare parameters
p_est = pp_params();
p_est = p_est.add_covar('rate',0,T_knots,'indicator'); % flat baseline rate
ps = cell(1,N_nodes);
ps{1} = p_est;
for n = 1:N_nodes-1
  ps{n+1} = ps{n}; ps{n+1}.response = n;
  ps{n+1} = ps{n}.add_covar('ens',n,R_knots,'spline');
end

% plot pair (later, sequence) of delay effect curves
for n = 1:N_nodes-1
  subplot(N_nodes,1,n+1);
  [t,y] = plot_spline(R_knots,R_b{n});
  plot(t,y); hold on;
end

% simulate data
% sim_dn = zeros(d.N_channels,d.T);
sim_dn = zeros(N_nodes,d.T);

% simulate poisson channel
sim_dn(1,:) = boolean(poissrnd(20*d.dt,[1,d.T]));
sd = pp_data(sim_dn,d.t);

% simulates other, delayed channels
m_true = pp_model();

for n = 2:N_nodes
  b = [T_b; reshape([R_b{1:n-1}],[],1)];
  m_true = m_true.make_X(sd,ps{n});
  sim_dn(n+1,:) = poissrnd(exp(m_true.X*b));
  sd.dn = sim_dn;
end
sd = pp_data(sim_dn,d.t);

% estimate delay effect curves
for n = 2:N_nodes
  m_est = pp_model();
  m_est = m_est.fit(sd,ps{n});
  R_b_est = m_est.b(ps{n}.covariate_ind{n});
  % compare estimate to ground truth
  [t,yhat] = plot_spline(R_knots,R_b_est);
  subplot(N_nodes-1,1,n-1), plot(t,yhat,'r');
end

%% rhythmic (AR) network

% set baseline rate
T_knots = [0 1];
T_b = log(20*d.dt);

% express delays as [knot,b] combos
Q_knots = [1 5 10 15];
Q_b = [-10 -3 0 1 0 0]';

% declare parameters
p_est = pp_params();
p_est.response = 1;
p_est = p_est.add_covar('rate',0,T_knots,'indicator'); % flat baseline rate
p_est = p_est.add_covar('self-hist',response,Q_knots,'spline'); % self-hist effects

% plot pair (later, sequence) of delay effect curves
[t,y] = plot_spline(Q_knots,Q_b);
figure, plot(t,y);
%%
% simulate data
% sim_dn = zeros(d.N_channels,d.T);
sim_dn = zeros(2,d.T);

% simulate poisson channel
sim_dn(2,:) = boolean(poissrnd(20*d.dt,[1,d.T]));
sd = pp_data(sim_dn,d.t);
burn_in = p_est.get_burn_in();

m_true = pp_model();
m_true = m_true.make_X(sd,p_est);

m_true.b = [T_b; Q_b];
N_cov = length(m_true.b);
X = zeros(1,N_cov);
X(1) = 1;

Xs = p_est.splineX(2);
sim_dn(1,1:burn_in) = poissrnd(exp(T_b),[1,burn_in]);

for t = burn_in+1:d.T
% % %   X = m_tru e.X(t,:); % this needs to be updated

  % use spike history to create X
  dn_T = sim_dn(response,t-[Q_knots(1):Q_knots(end)]);
  Xt = [1 dn_T*Xs];
  L = exp(X*m_true.b);
  sim_dn(1,t) = poissrnd(L);
end

sd = pp_data(sim_dn,d.t);
%%
% estimate self-history effect curves
m_est = pp_model();
m_est = m_est.fit(sd,p_est);
Q_b_est = m_est.b(p_est.covariate_ind{2});

% compare estimate to ground truth
[t,yhat] = plot_spline(Q_knots,Q_b_est);
hold on; plot(t,yhat,'r');

%% compute some cross-correlations for simulated data
figure
N=1;
N_lags=100;
% all_xc = zeros(N,2*N_lags+1);
for n = 1:N
  [xc,l] = xcorr(sd.dn(n,:),sd.dn(n+1,:),N_lags,'coeff');
%   all_xc(n,:) = xc;
  subplot(N,1,n), plot(l,xc);
end
%% compute auto-correlations
N=2;
N_lags=100;
% all_ac = zeros(N,N_lags+1);
for n = 1:N
  [ac,l] = autocorr(sd.dn(n,:)-mean(sd.dn(n,:)),N_lags);
%   all_ac(% figure, plot(l,all_xc(1,:))n,:) = ac;
  subplot(N,1,n), plot(l,ac); ylim([0,0.1]);
end