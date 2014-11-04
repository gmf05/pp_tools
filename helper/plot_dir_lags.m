function mov = plot_dir_lags(dir,lags)
  Fs = 3e4; % hard-coded sampling rate : change this if using for other purposes!
  Nl = size(dir,1);
  if nargin<2, lags = 1:Nl; end
  
  figure;
  mn = min(min(dir)); mx = max(max(dir));
  plot(dir(:,1),dir(:,2),'k');
  axis([mn,mx,mn,mx]);
  hold on;
  plot(0, 0, 'kx','markersize',2,'linewidth',2); hold on;
  for t = 1:Nl
    h1 = plot(dir(t,1),dir(t,2),'kx','markersize',10,'linewidth',2); 
    h2 = plot([0 dir(t,1)],[0 dir(t,2)],'b','linewidth',2); 
    if mod(t,20)==1, title([num2str(lags(t)) 'ms']); end
    pause(0.05); delete(h1,h2);
    mov(t) = getframe();
  end
end