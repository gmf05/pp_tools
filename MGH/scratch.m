% %%
count = 1;
% d = get_spikes('MG49','Seizure45','MUA');
d = get_spikes('MG49','Seizure45','LFP');
% d = get_spikes('MG49','Seizure45','ECoG');
p = pp_params();


% ms = cell(d.N_channels,6);
for response = 1:d.N_channels
% for response = 1
  p = pp_params();
  p = p.add_covar('rate',0,[0,1],'spline'); % baseline rate
%   p = p.add_covar('self-hist',response,[1 10:15:100 150:50:250 500],'spline');
%   p = p.add_covar('pop-hist',setdiff(1:d.N_channels,response),[0 10:15:100 150:50:250 500],'spline');
%   p.response = response;
%   p = p.add_covar('pop-up',c_up,[0 10:15:100 150:50:250 500],'spline');
%   p = p.add_covar('pop-down',c_down,[0 10:15:100 150:50:250 500],'spline');
%   p = p.add_covar('pop-left',c_left,[0 10:15:100 150:50:250 500],'spline');
%   p = p.add_covar('pop-right',c_right,[0 10:15:100 150:50:250 500],'spline');
  
  m = pp_model();
  m = m.fit(d,p); m.X=[];
  m
  ms{response,count} = m;
end
count = count+1;
