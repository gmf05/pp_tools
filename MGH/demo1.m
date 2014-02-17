%
load ../MG49_Seizure45_pp_thresh1_static0
cov_ind = 2; % self-history
% cov_ind = 3; % ensemble
cax = [0, 7];
Y = zeros(96,501+cov_ind-2);
Y0 = zeros(96,1);
for i = 1:96
  i
  if ~isempty(ms{i}.b)
    [t,y] = plot_spline(p.covariate_knots{cov_ind},ms{i}.b(p.covariate_ind{cov_ind}));
    Y(i,:) = exp(y); 
    Y0(i) = exp(y(1));
  end; 
end

Neuroport(Y,cax);
% pause;

%%
% cax = [0,10]
% Neuroport(Y0,cax);
