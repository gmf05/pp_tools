classdef Neuroport
  properties
    patient
    N_electrodes
    coord
    arrayMap
    Ws
    caxis
  end
  
  methods
    
    function obj = Neuroport(patient)
    obj.patient = patient;
    obj.N_electrodes = 96;
    
    % get coordinate remapping
    global DATA
    patient = obj.patient;
    num= xlsread([DATA '/' patient '/' patient '_lfp.xls']);
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
    
    function int = interior(obj)
      % Neuroport.interior:
      % Get interior electrodes
      % Usage: 
      % N = Neurport('MG49');
      % int = N.interior();
      %
      Xlo = min(obj.coord(:,1));
      Xhi = max(obj.coord(:,1));
      Ylo = min(obj.coord(:,2));
      Yhi = max(obj.coord(:,2));
      int = find(obj.coord(:,1)>Xlo & obj.coord(:,1)<Xhi & ...
                 obj.coord(:,2)>Ylo & obj.coord(:,2)<Yhi);
      % remove bad channels
      % for MG49: channel 89 is bad
      if isequal('MG49',obj.patient), int = setdiff(int, 89); end
    end
    
    function [cup,cdown,cleft,cright] = neighbors(obj,i)
      try
        cup = obj.arrayMap(obj.coord(i,1),obj.coord(i,2)+1);
        cdown = obj.arrayMap(obj.coord(i,1),obj.coord(i,2)-1);
        cleft = obj.arrayMap(obj.coord(i,1)-1,obj.coord(i,2));
        cright = obj.arrayMap(obj.coord(i,1)+1,obj.coord(i,2));
      catch 
        error('One (or more neighbors) do not exist');
      end      
    end
      
    function plot(obj,Ws,cax)    

    if nargin<2
      if isempty(obj.Ws)
        Ws = ones(96,1);
      else
        Ws = obj.Ws;
      end
    end
            
    if nargin<3 && isempty(obj.caxis)
      cax = [min(min(Ws)) max(max(Ws))];
    elseif nargin<3
      cax = obj.caxis;
    end    

    C = cax(2)-cax(1);    
    R = 0.5;
    if size(Ws,2)==obj.N_electrodes && size(Ws,1)~=obj.N_electrodes, Ws = Ws'; end
    T = size(Ws,2);

    % plotting---
    clf, set(gcf,'units','normalized','position',[0 0 1 1]);
    color_RGB = colormap();
%     my_img = phaseintensity(Ws);
    
    for t = 1:T
      for n = 1:obj.N_electrodes
        x = obj.coord(n,1);
        y = obj.coord(n,2);
        if C==0 || isnan(C)
          col_ind = 1;
        else          
          col_ind = round((Ws(n,t)-cax(1))/C*63)+1;
          col_ind = min(col_ind,64);
          col_ind = max(col_ind,1);
        end
        col = color_RGB(col_ind,:);
%         col = reshape(my_img(n,:,:),1,[]); % phase intensity
        fill([x-R x-R x+R x+R],[y-R y+R y+R y-R], col); hold on;
        text(x,y,num2str(n),'fontsize',22);
      end
      axis xy;
      set(gca,'XTick',[]);
      set(gca,'YTick',[]);
      
      if T>1
        title(num2str(t));
        pause(0.005);
        hold off;
      end
    end
  end
  
  function plot_dir(obj, bs, Ws)
    if nargin<3
      doCI = false;
    else
      doCI = true;
    end
    
    obj.plot();
    hold on;
    WT = [0 0 1 -1; -1 1 0 0];
    
    for n = 1:obj.N_electrodes
      if ~isempty(bs{response})
        b = bs{response};
        cumDir = WT*b;
        plot(obj.coord(n,1) + [-cumDir(1) 0], obj.coord(n,2) + [-cumDir(2) 0], 'b');
        if doCI
          cumDirW = WT * Ws{response} * WT';
          
        end
      end
    end
   end
   
  function [y0,W0] = spline_dir(obj, knots, b, W)
    % initialize variables
    s = 0.5;
%      WT = [0 0 1 -1; -1 1 0 0];
    Ns = length(knots) + 2;
    ys = cell(1,4);
    vr = cell(1,4);
    cv = cell(4,4); 
    
    % loop over the 4 directions, compute mean and variance
    % 1 = up->down (y-)
    % 2 = down->up (y+)
    % 3 = left->right (x+)
    % 4 = right->left (x-)
    for i = 1:4
      ind = (i-1)*Ns + (1:Ns);
      [lags, mn, mnp2se] = plot_spline(knots, b(ind), s, W(ind,ind));
      ys{i} = mn';
      vr{i} = ((mnp2se - mn)/2).^2;
      for j = 1:4, cv{i,j} = zeros(length(lags),1); end
    end
    
    N0 = length(lags);
    y0 = [ys{3} - ys{4}; ys{2} - ys{1}];
    W0 = zeros(2,2,N0);
    for n = 1:N0
      W0(1,1,n) = vr{3}(n) + vr{4}(n) - 2*cv{3,4}(n);
      W0(2,2,n) = vr{1}(n) + vr{2}(n) - 2*cv{1,2}(n);
      W0(1,2,n) = cv{1,2}(n) + cv{1,3}(n) + cv{1,4}(n) + cv{2,3}(n) + cv{2,4}(n) + cv{3,4}(n);
      W0(2,1,n) = W0(1,2,n);
    end
  end
  
  end
end
