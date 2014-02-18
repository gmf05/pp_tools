%% display P+I+E plots for all cells in a single session

date = '10-26-13';
protocol = 'prot0';
ind = 2;
s = Session(protocol,date,ind);
d = s.PPData;

% bad_ind = [1 6 10 15 25:27 31:32];
% good_ind = setdiff(1:d.N_channels, bad_ind);

figure('units','normalized','position',[0 0 1 1]);
% for response = good_ind
for response = 1:d.N_channels
  response
  load([d.Name '_c' num2str(response)]);
  m.plot(d,p);
  pause();
  for n = 1:3, subplot(3,1,n), hold on; end;
end
  

%% display P+I+E plots for a variety of cells
N_plots = 3;
protocol = 'prot1';
dispCells = {};

figure('units','normalized','position',[0 0 1 1]);
for i = 1:size(dispCells,1)
  date = dispCells{i,1};
  ind = dispCells{i,2};
  response = dispCells{i,3};
  s = Session(protocol,date,ind); d = s.PPData;
  filename = [d.Name '_c' num2str(response)]
  load(filename);
  m.plot(d,p);
  for n = 1:N_plots, subplot(N_plots,1,n), hold on; end;
end
  