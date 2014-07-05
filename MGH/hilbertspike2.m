function [spikeInd,amp] = hilbertspike2(D,thresh)

% get hilbert transform & its zero crossings (i.e. max/min of D)
H = imag(hilbert(D));
minInd = intersect(find(H(1:end-1)>0),find(H(2:end)<0));
maxInd = intersect(find(H(1:end-1)<0),find(H(2:end)>0));

% I = preliminary spike indices
I = [minInd maxInd];
Itype = [-1*ones(1,length(minInd)) ones(1,length(maxInd))];
[~,resort] = sort([minInd,maxInd]);
I = I(resort); Itype = Itype(resort);

% delete any max/max or min/min sequences
% (do these even occur??)
delInd = [];
for i = 2:size(I,2)-1
  if Itype(i-1)==Itype(i) || Itype(i+1)==Itype(i)
    delInd = [delInd i];
  end
end
I(delInd) = []; % drop bad sequences

% make sure we are starting with a max:
if I(1)==minInd(1), I(1) = []; end

% only keep indices where amplitude changes by 'thresh' units
zD = zscore(D); % assuming D is raw data here
spikeInd = [];
for i = 1:2:length(I)-1
  if zD(I(i))-zD(I(i+1))>=thresh
    [~,j] = min(zD(I(i):I(i+1)));
    spikeInd = [spikeInd I(i)+j-1];
%     minInd0 = [minInd0 I(i+1)];
  end
end

% STILL TO DO:
% eliminate "soft" minima that are found...


amp = D(spikeInd);

end