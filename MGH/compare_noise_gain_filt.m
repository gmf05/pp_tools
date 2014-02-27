% pp_tools;
% sz = seizure('MG49','Seizure45');
% d = sz.ECoG.PPData;
patient = 'MG49'; seizure = 'Seizure45'; 
load([DATA_DIR '/' patient '/' patient '_' seizure '_LFP_pp_thresh1']); d = obj; clear obj;
% d = sz.LFP.PPData;
% d = d.downsample(32);
p = pp_params();

response = 1;
p.response = response;

T_knots = [0 1];
p = p.add_covar('rate', 0, T_knots, 'indicator');

Q_knots = [1 21:20:101 151:50:501 1001];
p = p.add_covar('self-hist', response, Q_knots, 'spline');

R_knots = [1 21:20:101 151:50:501];
p = p.add_covar('ensemble', [1:response-1, response+1:d.N_channels], R_knots, 'spline');

% fit_method = 'glmfit'; noise = [];
fit_method = 'filt';
% fit_method = 'smooth';
% noise = [0 1e-4 1e-4]; % larger, ????
noise = [0 1e-6 1e-8]; % smaller, KS = 0.957, LL=-3.3720e03

all_KS = zeros(7,7,10);

for n1 = 1:7
	for n2 = 1:7

    noise(2) = 10^(-13+n1);
    noise(3) = 10^(-13+n2);

    p.fit_method = fit_method;
    p.noise = noise;
    p.downsample_est = 50;

    m = pp_model();

    % bs = cell(1, d.N_channels);
    for n = 1:10
	    response = 10*(n-1)+1;
	    p.response = response; p.covariate_channels{2} = response;
	    p.covariate_channels{3} = [1:response-1, response+1:d.N_channels];  
	    m = m.fit(d,p);
	    m
	    %   figure, m.plot(d,p); pause; hold off;
	    %   figure, m.gof(d); pause; close all;
	    %   bs{response} = m.b;
	    %   ms{response} = m;
      all_KS(n1,n2,n) = m.KS(1);
    end
  end
end

noise_range = -12:-6;
save compare_noise_gain_filt all_KS noise_range
