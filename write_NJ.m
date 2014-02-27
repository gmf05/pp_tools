switch result
	case 'interact'
    fpref = [date '-' protocol '-' num2str(ind) '_pp_interact_c'];
    I = 2;
  case 'spline'
    [];
  otherwise
    fpref = [date '-' protocol '-' num2str(ind) '_pp_c'];
    I = 1; 
end

% get file names
files = dir([fpref '*.mat']);
C = length(files);
noise_table = zeros(C,2);
laser_table = zeros(C,2);
interact_table = zeros(C,2);
KS_table = zeros(C,2);
max_chan = 0;

% loop over channels: load -> form table
for c = 1:C
  load(files(c).name);
	[i_on,i_off] = regexp(files(c).name, '(?<=_c)\d+');
	chan_str = files(c).name(i_on:i_off);
	chan = str2num(chan_str);
	max_chan = max(chan, max_chan);
  
  noise_table(chan,1) = m.b(end-I);
  noise_table(chan,2) = m.stats.p(end-I);
  laser_table(chan,1) = m.b(end-I+1);
  laser_table(chan,2) = m.stats.p(end-I+1);
%   interact_table(chan,1) = m.b(end-I+2);
%   interact_table(chan,2) = m.stats.p(end-I+2);
  KS_table(chan,:) = m.KS;
end

% show table
titles = {'channel', 'KS', 'KS (CI)', 'noise b', 'noise p', 'laser b', 'laser p'};
for i = 1:length(titles)-1
	fprintf([titles{i} '\t']);
end
fprintf([titles{end} '\n']);
table = [(1:max_chan)' KS_table noise_table laser_table];
table

% save table to csv file
csvFile = ['NJ_result_' result '.csv']; 
for i = 1:max_chan
	csvwrite(csvFile,date,i-1,0);
	csvwrite(csvFile,ind,i-1,1);	
end
csvwrite(csvFile,table,0,2);
