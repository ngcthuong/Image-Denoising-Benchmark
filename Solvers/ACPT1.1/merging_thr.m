function index_new=merging_thr(mean_cluster,indexes,k,thresh,X)
if k==1
    index_new=indexes;
else
    min2=zeros(1,k);
    min_i=zeros(1,k);
    for i=1:k
        ik=min(i,k-1);
        dif=sum((mean_cluster(:,ik)-mean_cluster(:,ik+1)).^2)/size(X,1);
        if ik==i
          min_i(i)=i+1;
        else
          min_i(i)=k-1; 
        end
        for j=i+1:k
                dif0 = sum((mean_cluster(:,i)-mean_cluster(:,j)).^2)/size(X,1);
                if dif0 ~= 0
                    if dif > dif0
                        dif = dif0;
                        min_i(i) = j;
                    end
                end
        end
        min2(1,i)=dif;
    end

    ind_re=1:k;
    for i = 1:k
        if min2(i)<thresh 
            ind_re(min_i(i))=min(i,ind_re(i));
        end
    end
    te=0;
    index_new={};
    for i=1:k
        II=find(ind_re==i);
        if size(II)
            te=te+1;
            in=[];
            for tt=1:length(II)
               in =[in ; indexes{II(tt)}];
            end
            index_new{te}=in;
        end
    end
end