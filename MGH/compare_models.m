%% fit model to big spikes
patient_name = 'MG49';
seizure_name = 'Seizure45';
data_type = 'LFP';
d1thresh = 0.5; d2thresh = 1;
lockout = 0.1; % time [sec]
N = Neuroport(patient_name);
int_elec = N.interior();
N_int = length(int_elec);

load MG49_window_all.mat
p = p_C;
T_knots = p.covariate_knots{1}; T_basis = p.covariate_bases{1}; T0 = length(T_knots);
Q_knots = p.covariate_knots{2}; Q_basis = p.covariate_bases{2}; Q = length(Q_knots);
R_knots = p.covariate_knots{3}; R_basis = p.covariate_bases{3}; R = length(R_knots);

% set params object 
d = get_spikes2(patient_name,seizure_name,data_type,d1thresh,d2thresh);
d = d.downsample(8); % do we want to downsample?

%%
j=2;
d0 = d.sub_time(tmins(j),tmaxs(j));
% for i=1:N_int
for i=1
  response=int_elec(i);
  p_N = change_response_MGH(p_N,response,N);
  p_S = change_response_MGH(p_S,response,N);
  p_I = change_response_MGH(p_I,response,N);
  p_C = change_response_MGH(p_C,response,N);
  m_N = ms_N{i,j};
  m_S = ms_S{i,j};
  m_I = ms_I{i,j};
  m_C = ms_C{i,j};
  DS = 0*d0.dn(response,:);
  [DN,LN]=pp_sim(m_N,d0,p_N,false); % 
  [DS,LS]=pp_sim(m_S,d0,p_S,false); % blew up with j=3, 1
  [DI,LI]=pp_sim(m_I,d0,p_I,true); % 
  [DC,LC]=pp_sim(m_C,d0,p_C,true);
end

D0 = pp_data([DN; DS; DI; DC; d0.dn(response,:)],d0.t);
figure, D0.plot('raster');


%%
figure(1), m_I.plot(d0,p_I);
figure(2), m_S.plot(d0,p_S);
figure(3), m_C.plot(d0,p_C);

%% hierarchical model comparison:

% compute delta Dev, Pearson's X^2, compare each with chi^2(p)
% ASSUME m1 = parent, m2 = child, m3 = child
% modelH = [0 1 1; 0 0 0; 0 0 0];
% [M,N] = find(modelH);

% for m = 1:length(M)
for i = 1:N_int
  for m = 1
%     parent = M(m); child = N(m);
%     mp = ms{parent};
%     mc = ms{child};
%     mp = ms_C{i,3};
%     mc = ms_I{i,3};
%     dDev = mc.dev - mp.dev;
%     dAIC = mc.AIC - mp.AIC;
%     P1 = size(mp.b,1); P2 = size(mc.b,1);
%     F1 = dDev/(P1-P2);
%     F2 = dAIC/(P1-P2);
%     pchi1 = 1-chi2cdf(dDev,P1-P2);
%     pchi2 = 1-chi2cdf(dAIC,P1-P2);
%     [int_elec(i), dDev, dAIC, pchi1,pchi2], pause;
%     pF1 = 1-fcdf(F1,P1,P1-P2);
%     pF2 = 1-fcdf(F2,P1,P1-P2);
    
    mC = ms_C{i,j}; mI = ms_I{i,j}; mS = ms_S{i,j};
    dDev_I = mI.dev - mC.dev;
    dDev_S = mS.dev - mC.dev;
    dAIC_I = mI.AIC - mC.AIC;
    dAIC_S = mS.AIC - mC.AIC;
    PC = size(mC.b,1); PI = size(mI.b,1); PS = size(mS.b,1);
    FCI = dDev_I/(PC-PI);
    FCS = dDev_S/(PC-PS);
    p_devCI = 1-chi2cdf(dDev_I,PC-PI);
    p_devCS = 1-chi2cdf(dDev_S,PC-PS);
    p_aicCI = 1-chi2cdf(dAIC_I,PC-PI);
    p_aicCS = 1-chi2cdf(dAIC_S,PC-PS);
    [int_elec(i),dDev_I,p_devCI,dDev_S,p_devCS,dAIC_I,p_aicCI,dAIC_S,p_aicCS]
    pause;
  end
end



%% test new version...
% [intrinsic, spatial, null, complete]
% modelH = [0 0 1 0; 0 0 1 0; 0 0 0 0; 1 1 1 0];
% modelH = [0 0 0 0; 1 0 0 0; 1 0 0 0; 1 1 1 0];
modelH = [0 0 1 0; 0 0 1 0; 0 0 0 0; 1 1 1 0];
% i=1; j=2;
for i=1:N_int
% ms = {ms_N{i,j} ms_I{i,j} ms_S{i,j} ms_C{i,j}};
ms = {ms_I{i,j} ms_S{i,j} ms_N{i,j} ms_C{i,j}};
ms{1}.name = 'intrinsic';
ms{2}.name = 'spatial';
ms{3}.name = 'null';
ms{4}.name = 'full';
table = compare_models2(modelH,ms)
pause;
end