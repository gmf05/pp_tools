PLOT_COLOR = 'b';
m = pp_model();
dt_ms = round(.001 / d.dt);

% T_knots, Q_knots, R_knots are assigned outside script
Q = length(Q_knots); R = length(R_knots);
chans = cellstr2num(d.labels);
int_elec = N.interior();
ps = [];
count = 1;
% for response = 1:d.N_channels
% for response = 21:40
% for response = 9
  response
  m = pp_model();
  p = pp_params();
  p.response = response;
  c_ind = chans(response);
  
  % first check that it's an interior electrode, then get up/down/etc
  if ismember(c_ind, int_elec)
    [c_up, c_down,c_left,c_right] = N.neighbors(c_ind);
    C_up = find(chans==c_up);
    C_down = find(chans==c_down);
    C_left = find(chans==c_left);
    C_right = find(chans==c_right);
    
    p = p.add_covar('rate',0,T_knots,T_basis); % baseline rate

    if Q>0
      p = p.add_covar('self-history',response,Q_knots,Q_basis);
    end
    
    if R>0
      p = p.add_covar('pop-hist1',C_up,R_knots,R_basis);
      if N_spatial_cov>1
        p = p.add_covar('pop-hist2',C_down,R_knots,R_basis);
        p = p.add_covar('pop-hist3',C_left,R_knots,R_basis);
        p = p.add_covar('pop-hist4',C_right,R_knots,R_basis);
      end
    end

    m = m.fit(d,p); m, %m.plot(d,p); pause(); %clf;
%     m = m.fit(d,p); m, m.gof(d); pause; clf;
    
    ps = p; % save parameters
  end
% end
p = ps;