% %%
% d = get_spikes('MG49','Seizure45','MUA');
d = get_spikes('MG49','Seizure45','LFP');
% d = get_spikes('MG49','Seizure45','ECoG');
% d = get_spikes('MG49','Seizure36','LFP');

ms = cell(d.N_channels,5);
ps = cell(d.N_channels,5);

count = 1;
for response = 1:d.N_channels
  p = pp_params();
  p.response = response;
  p = p.add_covar('rate',0,[0,1],'indicator'); % baseline rate
  ps{count} = p;
  
  m = pp_model();
  m = m.fit(d,p); m.X=[];
  m
  ms{response,count} = m;
end

count = count+1;
for response = 1:d.N_channels
  p = pp_params();
  p.response = response;
  
  p = p.add_covar('rate',0,[0,1],'indicator'); % baseline rate
  p = p.add_covar('self-hist',response,[1 10:15:100 150:50:250 500],'spline');
  ps{count} = p;
  
  m = pp_model();
  m = m.fit(d,p); m.X=[];
  m
  ms{response,count} = m;
end

count = count+1;
for response = 1:d.N_channels
  p = pp_params();
  p.response = response;
  
  p = p.add_covar('rate',0,[0,1],'indicator'); % baseline rate
  p = p.add_covar('pop-hist',setdiff(1:d.N_channels,response),[0 10:15:100 150:50:250 500],'spline');
  ps{count} = p;
  
  m = pp_model();
  m = m.fit(d,p); m.X=[];
  m
  ms{response,count} = m;
end

count = count+1;
for response = 1:d.N_channels
  p = pp_params();
  p.response = response;
  
  p = p.add_covar('rate',0,[0,1],'indicator'); % baseline rate
  p = p.add_covar('self-hist',response,[1 10:15:100 150:50:250 500],'spline');
  p = p.add_covar('pop-hist',setdiff(1:d.N_channels,response),[0 10:15:100 150:50:250 500],'spline');
  ps{count} = p;
  
  m = pp_model();
  m = m.fit(d,p); m.X=[];
  m
  ms{response,count} = m;
end

%%
% load mg49_s45_lfp_models
N = Neuroport('MG49');
count = 5;
for response = 1:d.N_channels
  p = pp_params();
  p.response = response;
  c_ind = str2num(d.Labels{response});
  ms{response,count} = [];
  % first check that it's an interior electrode, then get up/down/etc
  if max(N.coord(c_ind,:))<10 && min(N.coord(c_ind,:))>1
    % match Labels{response} to c_ind
    c_ind = str2num(d.Labels{response});
    c_up = N.arrayMap(N.coord(c_ind,1),N.coord(c_ind,2)+1);
    c_down = N.arrayMap(N.coord(c_ind,1),N.coord(c_ind,2)-1);
    c_left = N.arrayMap(N.coord(c_ind,1)-1,N.coord(c_ind,2));
    c_right = N.arrayMap(N.coord(c_ind,1)+1,N.coord(c_ind,2));
    % get indices in labels
    
    clear C_up C_down C_left C_right
    for j = 1:d.N_channels
      if str2num(d.Labels{j})==c_up, C_up = j; end
      if str2num(d.Labels{j})==c_down, C_down = j; end
      if str2num(d.Labels{j})==c_left, C_left = j; end
      if str2num(d.Labels{j})==c_right, C_right = j; end
    end
    try
    p = p.add_covar('rate',0,[0,1],'indicator'); % baseline rate
    p = p.add_covar('self-hist',response,[1 10:15:100 150:50:250 500],'spline');
    p = p.add_covar('pop-hist1',C_up,[0 20 100 150 500],'spline');
    p = p.add_covar('pop-hist2',C_down,[0 20 100 150 500],'spline');
    p = p.add_covar('pop-hist3',C_left,[0 20 100 150 500],'spline');
    p = p.add_covar('pop-hist4',C_right,[0 20 100 150 500],'spline');

    m = pp_model();
    m = m.fit(d,p); m.X=[];
    m
    [c_ind,c_up,c_down,c_left,c_right]
%     ms{count} = m;
    ms{response,count} = m;
    end
  end
end
%%

save mg49_s45_lfp_models ms ps
% save mg49_s36_lfp_models ms p
