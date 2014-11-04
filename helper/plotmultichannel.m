function plotmultichannel(t,d)

for i = 1:size(d,1)
  plot(t,normalize(zscore(d(i,:)))-0.5+i); hold on
end

end