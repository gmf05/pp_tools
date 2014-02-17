function Neuroport(Ws,cax)
% Y = zeros(96,501);
% for i = 1:96, i, if ~isempty(ms{i}.b), [t,y] = plot_spline(p.covariate_knots{2},ms{i}.b(p.covariate_ind{2})); Y(i,:) = exp(y); end; end

if nargin<1, Ws = ones(96,1); end
if nargin<2
  if max(range(Ws))==0
    cax = [Ws(1)-1 Ws(1)];
  else
    cax = [min(min(Ws)) max(max(Ws))];
  end
end

C = cax(2)-cax(1)
R = 0.5;
N_electrodes = 96;
[T,temp] = setdiff(size(Ws),N_electrodes);
% make sure columns of Ws are spatial weights:
if temp==1, Ws = Ws'; end; clear temp;

% set coordinates
coord = zeros(N_electrodes,2);
count = 1;
for i = 1:10
  if i==1 || i==10
    for j = 2:9
      coord(count,:) = [j,i];
      count = count+1;
    end
  else
    for j = 1:10
      coord(count,:) = [j,i];
      count = count+1;
    end
  end
end

% plotting---

% initialize figure
figure('units','normalized','position',[0 0 1 1]);

% set up colormap
colormap('default');
color_RGB = colormap();

for t = 1:T
  for n = 1:N_electrodes
    % set coordinates
    x = coord(n,1);
    y = coord(n,2);
    % get color
    col_ind = round((Ws(n,t)-cax(1))/C*63)+1;
    col_ind = min(col_ind,64); col_ind = max(col_ind,1);
    col = color_RGB(col_ind,:);
    % fill    
    fill([x-R x-R x+R x+R],[y-R y+R y+R y-R], col); hold on;
    text(x,y,num2str(n),'fontsize',20); hold on;
  end
  axis ij;
%   title(num2str(t*0.011*1e3))
  title(num2str(t))
  pause(0.6);
  hold off;
end

set(gca,'XTickLabel',{});
set(gca,'YTickLabel',{});


end
