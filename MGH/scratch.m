%%

diary small_spike_models.txt

DO_CONF_INT = false;
T_knots = p.covariate_knots{1}; T_basis = p.covariate_bases{1};
Q_knots = p.covariate_knots{2}; Q_basis = p.covariate_bases{2}; 
R_knots = p.covariate_knots{3}; R_basis = p.covariate_bases{3};
Q = length(Q_knots); R = length(R_knots);

m0 = pp_model();
% m0.b = b;
% m0.W = W;
m0.fit_method = 'glmfit';
chans = cellstr2num(d.labels);
int_elec = N.interior();

for response = 1:d.N_channels
  
  response
  p0 = pp_params();
  p0.response = response;
  c_ind = chans(response);
  
  if ismember(c_ind,int_elec)
    [c_up, c_down,c_left,c_right] = N.neighbors(c_ind);
    C_up = find(chans==c_up);
    C_down = find(chans==c_down);
    C_left = find(chans==c_left);
    C_right = find(chans==c_right);
  %     D0 = d.sub_data([response,C_up,C_down,C_left,C_right]);

    p0 = p0.add_covar('rate',0,T_knots,T_basis); % baseline rate

    if Q>0
      p0 = p0.add_covar('self-history',response,Q_knots,Q_basis);
%       p0 = p0.add_covar('self-history2',response+1,Q_knots,Q_basis);
    end

    if R>0
      p0 = p0.add_covar('pop-hist1',C_up,R_knots,R_basis);
      p0 = p0.add_covar('pop-hist2',C_down,R_knots,R_basis);
      p0 = p0.add_covar('pop-hist3',C_left,R_knots,R_basis);
      p0 = p0.add_covar('pop-hist4',C_right,R_knots,R_basis);
    end

    % m0 = m0.makeX(d,p0);
    m0 = m0.fit(d,p0);
    if m0.KS(1)>0.1, PLOT_COLOR = 'r';
    else PLOT_COLOR = 'b'; end
    m0.plot(d,p0); hold on; pause(0.2);
    m0
  end
  
end

% m0.y = d.dn(response,:)';
% m0.CIF = exp(m0.X*b);
% m0 = m0.calcGOF();
diary off;
