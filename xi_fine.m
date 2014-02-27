% pt_name = 'MG49';  sz_name = 'Seizure45';
% 
% load([pt_name '_' sz_name '_pp_thresh1_static.mat']);
N_electrodes = 96;
ms2 = cell(1,N_electrodes);
% d must be high-res data
% sz = seizure(pt_name,sz_name);
% d = sz.LFP.PPData;

for i = 1:N_electrodes
  i
  p.response = i;
  p.covariate_channels{2} = i;
  p.covariate_channels{3} = [1:i-1 i+1:N_electrodes];
  if ~isempty(ms{i}.b)
    m = ms{i}.make_X(d,p); X = m.X; clear m;
    ms2{i} = pp_model();
    ms2{i}.CIF = glmval(ms{i}.b,X,'log','constant','off');
  end
end

XI_fine = xi(ms2,d);
save([pt_name '_' sz_name '_groups_XI_fine'],'XI_fine','p');
