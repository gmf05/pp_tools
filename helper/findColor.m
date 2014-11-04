function rgb = findColor(X,cax)

if nargin<2, cax = [min(min(X)) max(max(X))]; end
minX = cax(1); 
rangeX = range(cax);

N = length(X);
rgb = zeros(N,3);
cmap = colormap();

for n = 1:N
  i = round((X(n)-minX)/rangeX*63)+1;
  rgb(n,:) = cmap(i,:);
end

end