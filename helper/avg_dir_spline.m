function [Ds, Derr] = avg_dir_spline(b,cov_knots,spatial_cov_ind)
Nl = range(cov_knots)+1;
Xs = zeros(Nl,1);
Ys = zeros(Nl,1);
dir = [0 -1; 0 1; 1 0; -1 0];
for i = 1:4
  [~,y] = plot_spline(cov_knots,b(spatial_cov_ind(:,i)));
  Xs = Xs + y.*dir(i,1);
  Ys = Ys + y.*dir(i,2);
end
Ds = [Xs Ys];
Derr = [];
end