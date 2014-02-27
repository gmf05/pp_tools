load MG49_model_hierarchy;
% load hi-res data
load /projectnb/ecog/Data/MG49/MG49_Seizure36_LFP_pp_thresh1_hires
d = data; clear data;

N_channels = 96;
Xis = zeros(N_channels,N_channels,4);
mnull = cell(1,N_channels);
mens0 = cell(1,N_channels);
mens1 = cell(1,N_channels);
mens2 = cell(1,N_channels);

ens1 = [1:36 38 42 46:2:52 58:61 63 95];
ens2 = setdiff(1:N_channels,ens1);

m = pp_model();
fprintf(['Evaluating conditional intensity functions...\n']);
for i = 1:N_channels
  i
  % assign mnull, mens0 etc
  p_null.response = i;
  p_ens0.response = i;
  p_ens1.response = i;
  p_ens2.response = i;

  p_ens0.covariate_channels{2} = i;
  p_ens1.covariate_channels{2} = i;
  p_ens2.covariate_channels{2} = i;

  p_ens1.covariate_channels{3} = setdiff(1:N_channels,i);
  p_ens2.covariate_channels{4} = setdiff(ens1,i);
  p_ens2.covariate_channels{4} = setdiff(ens2,i);

  m = pp_model();
  mnull{i} = m;
  mens0{i} = m;
  mens1{i} = m;
  mens2{i} = m;

  if ~isempty(ms_ens2{i,1})
    m = m.make_X(d,p_null);
    mnull{i}.CIF = glmval(ms_null{i,1}.b,m.X,'log','constant','off'); m.X=[];
    m = m.make_X(d,p_ens0);
    mens0{i}.CIF = glmval(ms_ens0{i,1}.b,m.X,'log','constant','off'); m.X=[];
    m = m.make_X(d,p_ens1);
    mens1{i}.CIF = glmval(ms_ens1{i,1}.b,m.X,'log','constant','off'); m.X=[];
    m = m.make_X(d,p_ens2);
    mens2{i}.CIF = glmval(ms_ens2{i,1}.b,m.X,'log','constant','off'); m.X=[];
  end
end

fprintf('Computing xi...\n');
Xis(:,:,1) = xi(mnull,d);
fprintf('Null done\n');
Xis(:,:,2) = xi(mens0,d);
fprintf('Ens0 done\n');


save MG49_Seizure36_Xi_hierarchy Xis

Xis(:,:,3)= xi(mens1,d);
fprintf('Ens1 done\n');
Xis(:,:,4) = xi(mens2,d);
fprintf('Ens2 done\n');

save MG49_Seizure36_Xi_hierarchy Xis
