dates = {'5-28-13', '6-11-13', '6-13-13', '6-14-13', '6-18-13', '7-12-13', '7-18-13', '7-23-13', '7-24-13', '10-3-13', '10-10-13', '10-17-13', '12-22-13'};
protocol = 'prot1';

% for d_ind = 1:length(dates)
% for d_ind = 10:length(dates)
for d_ind = 6:9
	date = dates{d_ind};
	for ind = 1:2
		run_NJ_interact;
	end
end
