function XI = xi(ms,d,lag)
  if nargin<3, lag = 0; end

  cifs = zeros(d.N_channels,length(ms{1}.CIF));
  for i = 1:length(ms)
    if ~isempty(ms{i}) && ~isempty(ms{i}.CIF)
      cifs(i,:) = ms{i}.CIF;
    end
  end
  
  XI = zeros(d.N_channels);
  
  if lag==0
    for i = 1:d.N_channels
      i
      for j = i+1:d.N_channels
        numer = sum(d.dn(i,:).*d.dn(j,:));
        denom = sum(cifs(i,:).*cifs(j,:));
        if denom>0
          XI(i,j) = numer / denom;
        end
      end
    end
    XI = XI + XI';
  
  else
    for i = 1:d.N_channels
      i
      for j = 1:d.N_channels
        numer = sum(d.dn(i,lag+1:end).*d.dn(j,1:end-lag));
        denom = sum(cifs(i,lag+1:end).*cifs(j,1:end-lag));
        if denom>0
        	XI(i,j) = numer / denom; 
        end
      end
    end
  end
  
  XI(XI==0) = 1;
end
