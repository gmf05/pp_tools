function rgb = interpolateColor(c1,c2,steps)
if nargin<3, steps = 64; end
rgb = zeros(steps,3);
for n = 1:3
  rgb(:,n) = c1(n) + (c2(n)-c1(n))/(steps-1) * (0:steps-1);
end
end