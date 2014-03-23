classdef Neuroport
  properties
    patient
    N_electrodes
    coord
    arrayMap
    Ws
    color
  end
  
  methods
    
    function obj = Neuroport(patient)
    obj.patient = patient;
    obj.N_electrodes = 96;
    
    % get coordinate remapping
    global DATA_DIR
    patient = obj.patient;
    num= xlsread([DATA_DIR '/' patient '/' patient '.xls']);
    arrayMap = NaN(10,10);
    for n = 1:96;
      arrayMap(num(n,1)+1,num(n,2)+1) = n;
    end
    obj.arrayMap = arrayMap;
    
    % set coordinates
    coord = zeros(obj.N_electrodes,2);
    for i = 1:10
      if i==1 || i==10
        col_ind = 2:9;
      else
        col_ind = 1:10;
      end

      for j = col_ind
        coord(arrayMap(j,i),:) = [j,i];
      end
    end
    obj.coord = coord;
  end
  
  function plot(obj,Ws,cax)
    
    % Y = zeros(96,501);
    % for i = 1:96, i, if ~isempty(ms{i}.b), [t,y] = plot_spline(p.covariate_knots{2},ms{i}.b(p.covariate_ind{2})); Y(i,:) = exp(y); end; end

    if nargin<2
      if isempty(obj.Ws)
        Ws = ones(96,1);
      else
        Ws = obj.Ws;
      end
    end
    
    if isempty(obj.color) && nargin<3
      if max(range(Ws))==0
        cax = [Ws(1)-1 Ws(1)];
      else
        cax = [min(min(Ws)) max(max(Ws))];
      end
    elseif isempty(cax)
      cax = obj.color;
    end

    C = cax(2)-cax(1);
    R = 0.5;
    if size(Ws,2)==obj.N_electrodes && size(Ws,1)~=obj.N_electrodes, Ws = Ws'; end
    T = size(Ws,2);

    % plotting---
    figure('units','normalized','position',[0 0 1 1]);
    
    % set up colormap
    colormap('default'), color_RGB = colormap();

    for t = 1:T
      for n = 1:obj.N_electrodes
        % set coordinates
        x = obj.coord(n,1);
        y = obj.coord(n,2);
        % get color
        col_ind = round((Ws(n,t)-cax(1))/C*63)+1;
        col_ind = min(col_ind,64); col_ind = max(col_ind,1);
        col = color_RGB(col_ind,:);
        % fill
        fill([x-R x-R x+R x+R],[y-R y+R y+R y-R], col); hold on;
        text(x,y,num2str(n),'fontsize',22);
      end
      axis xy;
      title(num2str(t));
      pause(0.005);
      hold off;
      set(gca,'XTick',[]);
      set(gca,'YTick',[]);
    end
  end
  end
end