%%

for i = 25:96
  i
  fname = [DATA_DIR '/MG49/MG49_filt/MG49_Seizure45_LFP_pp_filt_c' num2str(i)];
  open([fname '.fig']);
  subplot(311); title(['Electrode ' num2str(i)]);
  subplot(312); ylim([0,500]); caxis([0,2]); colorbar;
  subplot(313); ylim([0,500]); caxis([0.9,1.1]); colorbar;
  pause();
%   plot2svg([fname '.svg']);
  close;
end


%%
% s = seizure('BU1','Seizure1');
s = seizure('MG49','Seizure45');
d = s.ECoG.PPData;
p = pp_params();

response = 1;
p.response = response;

T_knots = [0 1];
p = p.add_covar('rate', 0, T_knots, 'indicator');

Q_knots = [1 21:20:101 151:50:501 1001];
p = p.add_covar('self-hist', response, Q_knots, 'spline');

R_knots = [1 21:20:101 151:50:501];
p = p.add_covar('ensemble', [1:response-1, response+1:d.N_channels], R_knots, 'spline');

fit_method = 'glmfit'; noise = [];
% fit_method = 'filt';
% fit_method = 'smooth';
% noise = [0 1e-8 1e-10]; % small, seems to work well
noise = [0 1e-6 1e-8]; % small, seems to work well
% noise = [0 1e-8 1e-10];
% noise = [0 1e-8 1e-10]; 

p.fit_method = fit_method;
p.noise = noise;
p.downsample_est = 5;

m = pp_model();
m = m.fit(d,p);
m

% figure, m.plot(d,p);
% figure, m.gof(d);

