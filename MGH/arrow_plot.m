function [x0,y0] = arrow_plot(m,p,ctr)

global PLOT_COLOR

% translate center, if desired
if nargin<3, ctr = [0 0]; end

% get [up,down,left,right] splines
Ncov = p.covariate_ind{1}(end);
Nlags = p.covariate_knots{1}(end)-p.covariate_knots{1}(1)+1;
Y = zeros(4,Nlags);
for n = 1:4
  [t,y] = plot_spline(p.covariate_knots{1},m.b(Ncov*(n-1)+(1:Ncov)));
  Y(n,:) = exp(y);
end

% direction is a sum over splines:
% (result is 2 vectors over lag times)
y0 = Y(2,:)-Y(1,:);
x0 = Y(3,:)-Y(4,:);

% % % plot over lag times
% % plot(x0,y0); hold on;
% % for i = 1:length(x0)
% %   h=plot(ctr(1)+[0 x0(i)],ctr(2)+[0 y0(i)],'r');
% %   pause(0.05);
% %   delete(h);
% % end

% plot single arrow
% summarize lags into single number...?
%
% then plot...
SCALE = 0.1;
x = x0(1)*SCALE;
y = y0(1)*SCALE;

plot(ctr(1)+[0 x], ctr(2)+[0 y], PLOT_COLOR);

end