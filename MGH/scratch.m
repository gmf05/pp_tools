
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

%% compute deviance & AIC for models missing these measures

for i = 1:96
  i
  if ~isempty(ms_ens2{i,1})
    ms_null{i,1}.dev = 2*(sum(log(poisspdf(ms_null{i,1}.y,ms_null{i,1}.y))) - ms_null{i,1}.LL);
    ms_null{i,1}.AIC = ms_null{i,1}.dev + 2*length(ms_null{i,1}.b);
    
    ms_ens0{i,1}.dev = 2*(sum(log(poisspdf(ms_ens0{i,1}.y,ms_ens0{i,1}.y))) - ms_ens0{i,1}.LL);
    ms_ens0{i,1}.AIC = ms_ens0{i,1}.dev + 2*length(ms_ens0{i,1}.b);

    ms_ens1{i,1}.dev = 2*(sum(log(poisspdf(ms_ens1{i,1}.y,ms_ens1{i,1}.y))) - ms_ens1{i,1}.LL);
    ms_ens1{i,1}.AIC = ms_ens1{i,1}.dev + 2*length(ms_ens1{i,1}.b);
    
    ms_ens2{i,1}.dev = 2*(sum(log(poisspdf(ms_ens2{i,1}.y,ms_ens2{i,1}.y))) - ms_ens2{i,1}.LL);
    ms_ens2{i,1}.AIC = ms_ens2{i,1}.dev + 2*length(ms_ens2{i,1}.b);
  end
  
  if ~isempty(ms_ens2{i,2})
    ms_null{i,2}.dev = 2*(sum(log(poisspdf(ms_null{i,2}.y,ms_null{i,2}.y))) - ms_null{i,2}.LL);
    ms_null{i,2}.AIC = ms_null{i,2}.dev + 2*length(ms_null{i,2}.b);
    
    ms_ens0{i,2}.dev = 2*(sum(log(poisspdf(ms_ens0{i,2}.y,ms_ens0{i,2}.y))) - ms_ens0{i,2}.LL);
    ms_ens0{i,2}.AIC = ms_ens0{i,2}.dev + 2*length(ms_ens0{i,2}.b);
    
    ms_ens1{i,2}.dev = 2*(sum(log(poisspdf(ms_ens1{i,2}.y,ms_ens1{i,2}.y))) - ms_ens1{i,2}.LL);
    ms_ens1{i,2}.AIC = ms_ens1{i,2}.dev + 2*length(ms_ens1{i,2}.b);
    
    ms_ens2{i,2}.dev = 2*(sum(log(poisspdf(ms_ens2{i,2}.y,ms_ens2{i,2}.y))) - ms_ens2{i,2}.LL);
    ms_ens2{i,2}.AIC = ms_ens2{i,2}.dev + 2*length(ms_ens2{i,2}.b);
  end
end