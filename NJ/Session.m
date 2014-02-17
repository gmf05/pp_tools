classdef Session
  properties
    Name
    noise_TS
    laser_TS    
    PPData
  end
  
  methods
    % Constructor    
    function obj = Session(protocol_name, date, ind)
      
      fprintf(['\n\nLoading point process data from session ' ...
        num2str(ind) ' on ' date ', ' protocol_name '\n']);      
      
      global DATA_DIR
      OLD_DIR = pwd();
      cd([DATA_DIR '/NJ/' protocol_name]);
      raw_filename = [date '-' protocol_name '-' num2str(ind)];
      pp_filename = [raw_filename '_pp'];
      
       if exist([pp_filename '.mat'],'file')
        fprintf(['Point process data found for session ' ...
          num2str(ind) ' on ' date '.\n']);
        load(pp_filename);
       else
        fprintf(['Point process data not found for session ' ...
             num2str(ind) ' on ' date  '.\nCreating it...\n']);
        % create pp_data
        load(raw_filename);
        
        % incorporate laser TS, noise TS????        
        tmin = spike_data.tbeg;
        tmax = spike_data.tend;
%         tmax = 300;
        dt = 1/spike_data.freq;
        t_axis = tmin:dt:tmax;
        T = length(t_axis);
        
        N_channels = length(spike_data.neurons);
        labels = cell(1,N_channels);
        dn = zeros(N_channels, T);
        
        for n = 1:N_channels
          spikes = spike_data.neurons{n}.timestamps;
          spikes = spikes(spikes>=tmin & spikes<=tmax);
          dn(n,:) = hist(spikes,t_axis);
          labels{n} = spike_data.neurons{n}.name;
        end
        d = pp_data(dn,t_axis,pp_filename,labels);
        
        obj.Name = pp_filename;
        obj.PPData = d;
        obj.noise_TS = noise_TS;
        obj.laser_TS = laser_TS;
        save([pp_filename '.mat'],'obj','-v7.3');
      end
      
      
      cd(OLD_DIR);
      fprintf('Done!\n\n');
    end
      
  end
end
      
  
  
