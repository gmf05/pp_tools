
Ne = 96; % number of electrodes
Ws = zeros(Ne,1);
cov_ind = 2;
bad_val = 0;
Y = zeros(Ne,501);
for i = 1:Ne
  try
    [t,y] = plot_spline(p.covariate_knots{cov_ind},ms{i}.b(p.covariate_ind{cov_ind}));
%     Y(i,:) = exp(y);
%     Ws(i) = exp(max(y));
%     Ws(i) = exp(y(1));
%     Ws(i) = ms{i}.KS(1);
  catch    
%     Ws(i) = bad_val;
  end
end

%%

% cax = [min(Ws),max(Ws)]; % color axis
Neuroport(Ws,cax);
