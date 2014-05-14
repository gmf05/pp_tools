classdef patient
  properties
    Name
    Comments
    Seizures
    Sex
    AgeOnset
    AgeSurgery
    NumSeizures
    Etiology
    Lobe
    Engel
  end
  
  methods
    % Constructor
    function obj = patient(Name); 
      if nargin<1, Name = '';
      else
        global DATA_DIR
        OLD_DIR = pwd();
        cd([DATA_DIR '/' Name]);
        eval(['obj = ' Name '();']);
        cd(OLD_DIR);
      end
    end
    
    function d = dir(obj)
      global DATA_DIR
      d = [DATA_DIR '/' obj.Name];
    end
    
    function GUI(obj)
      global GUI_FIG;
      GUI_FIG = figure();      
      
      % draw LFP, ECoG click buttons
      % 
      
      subplot('position',[0.1,0.1,0.8,0.7]);
      obj.plot('ECoG');
    end
    
    function plot(obj,data_type,Ws,W_axis)
      if nargin<3
        Weights = false;
        Ws = [];
        W_axis = [];        
      else
        Weights = true;
        N_electrodes = length(Ws);
        if nargin<4, W_axis = [min(Ws) max(Ws)]; end
      end
      
      % set up colormap for plotting Ws      
      if Weights
        colormap('default');
        caxis(W_axis);
        color_RGB = colormap();
      else
        global PLOT_COLOR % if weights are uniform
        col = PLOT_COLOR;
      end
      
      global ELEC_COORD;
      switch data_type
        case 'ECoG'
          patient_img = imread([obj.dir() '/' obj.Name '.png']);
          imagesc(patient_img), hold on;
          colormap(gray); freezeColors;
          load([obj.dir() '/' obj.Name '_ECoG_map']);
          ELEC_COORD = ECoG.Coord;
        case 'LFP'
          load([obj.dir() '/' obj.Name '_LFP_map']);
          ELEC_COORD = LFP.Coord;
      end
      N_electrodes = length(ELEC_COORD);
      R = 8;
      
      for n = 1:N_electrodes
        
        if Weights
          % SOMETHING SEEMS WRONG WITH W_ind -- too small?
          [~,W_ind] = min(abs(W_axis - Ws(n)));
          col = color_RGB(W_ind,:);
        end
        
        x = ELEC_COORD(n,1);
        y = ELEC_COORD(n,2);
        xs = [x-R/2 x-R/2 x+R/2 x+R/2];
        ys = [y-R/2 y+R/2 y+R/2 y-R/2];
        fill(xs,ys,col); hold on;
        text(10*n,10*n,num2str(n),'k','LineWidth',10,'fontsize',100);
      end
      
      set(gca,'YTick',[]); 
      set(gca,'XTick',[]);
      
%       update_fig();
      
    end
  end
end

      

      
      