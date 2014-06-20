function [multiW, multiY, multiL] = spline_cov(multiknots,b,W,s)
% % I = p.covariate_ind{2}(1):p.covariate_ind{end}(end); W0 = spline_cov({p.covariate_knots{2:end}},m.b(I),m.W(I,I));

    % checks::
    if nargin<4 || isempty(s), s=0.5; end;
    
        % tension matrix
    s_coeff = [-s  2-s s-2  s; 2*s s-3 3-2*s -s; ...
               -s   0   s   0;   0   1   0   0];
    
    Nk = length(multiknots);    
    % initialize multiX
    multiX = [];
    multiY = cell(1,Nk);
    multiL = cell(1,Nk);
    ind = 0;
    for k = 1:Nk
      knots = multiknots{k};
      ind = ind(end) + (1:length(knots)+2);
      if  range(knots)>1000, dtau = 2;
      elseif range(knots)>1, dtau = 1; 
      elseif range(knots) dtau = 1e-2; end;
      
      % covariate axis
      tau = knots(1):dtau:knots(end);
      multiL{k} = tau;
      NT = length(tau);

      % set X0 based on knots        
      X0 = zeros(NT,length(knots)+2);    
      spacing = diff(knots)*1/dtau;
      count=1;
      for i=1:length(spacing)
          alphas = (0:spacing(i)-1)./spacing(i);
          X0(count:count+spacing(i)-1, i:i+3) = [alphas'.^3 alphas'.^2 alphas' ones(spacing(i),1)] * s_coeff;
          count = count+spacing(i);
      end
      X0(end, i:i+3) = [1 1 1 1] * s_coeff; % alpha = 1
      Nr = size(X0,1); Nc = size(X0,2);      
      multiX(end+(1:Nr),end+(1:Nc)) = X0;
      multiY{k} = (X0*b(ind))';
    end
    
    multiW = multiX*W*multiX';
end