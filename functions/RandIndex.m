function [RI] = RandIndex(Y_star,Y_true)

    [mYs,nYs]=size(Y_star);
    [mYt,nYt]=size(Y_true);
    TP=0;TN=0;FN=0;FP=0;
    for i=1:nYs
        for j=i+1:nYs
            n1=sum(Y_star(:,i)==Y_star(:,j));
            n2=sum(Y_true(:,i)==Y_true(:,j));
            if ((n1==mYs)&&(n2==mYt))
                TP=TP+1;
            elseif((n1<mYs)&&(n2<mYt))
                TN=TN+1;
            elseif((n1==mYs)&&(n2<mYt))
                FN=FN+1;
            else
                FP=FP+1;
            end
        end
    end
    
    RI = (TP+TN)/(TP+TN+FN+FP);
    end
