function inds=getclosest(data,vals)
        inds=zeros(1,length(vals));
        
        for i=1:length(vals)
            [y ind]=min(abs(data-vals(i)));
            inds(i)=ind;
        end

end