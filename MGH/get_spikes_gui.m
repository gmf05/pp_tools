function spikeInd = get_spikes_gui(data,t,tstart,tstop)

dt = t(2)-t(1);
plot(t,data); hold on;
xlim([tstart,tstop]);

% spikeTimes = [];
spikeInd = [];
doSpikeFind = true;

while doSpikeFind

%   disp('(S)pike find / (L)eft / (R)ight? ');
  query1 = getkey();
  switch query1
    case 28 % left arrow
      scroll(-1,1); pause(0.1);
    case 29 % right arrow
      scroll(1,1); pause(0.1);
    case 115 % s
      % get user clicks
      % turn into interval bounds
      A = getrect(); % xmin ymin width height
      ion = getclosest(t,A(1));
      ioff = ion + round(A(3)/dt);

      % find min on interval
      [~,minind] = min(data(ion:ioff));
      tind = ion + minind - 1;

      % display x on putative spike
      h = plot(t(tind), data(tind),'rx','markersize',8,'linewidth',3);

      % ask for user confirmation
      disp('Is this a true spike? (y/n): ');
      query2 = getkey();
      switch query2
        case 121 % y
          spikeInd = [spikeInd tind];        
%         case 110 % n
        otherwise
          delete(h);
      end
    case 113 % q
      doSpikeFind = false;
    otherwise
      pause();
  end
end

% spikeTimes = sort(spikeTimes);
spikeInd = sort(spikeInd);

end