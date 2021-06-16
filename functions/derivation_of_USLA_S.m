function Deri_of_S=derivation_of_USLA_S(X,ST_t,Shape_t,Xkj_skl,SSij_sil,SS,Parameter)
DShape_t=Shape_t(:,2:end);
[mST,nST]=size(ST_t);
[mDShape,nDShape]=size(DShape_t);
%Der_of_S=zeros(mDShape,nDShape);

%Part1 YLY'
%Parameter.sigma1 = 1;

parameter2=Parameter.lambda_1*SS;
for k=1:mDShape
    len=Shape_t(k,1);
    for l=1:len
        for i=1:mST
            for j=1:nST
                P1(i,j)=ST_t(i,j)*(1/Parameter.sigma1^2)*(X(k,i)-X(k,j))*(Xkj_skl(k,i,l)-Xkj_skl(k,j,l));
%                 P1(i,j)=STij_skl(i,j,k,l);
            end
%             P1(i,j)=STij_skl(i,j,k,l);
        end
       Part1(k,l)=sum(sum(P1));
       %Part2 
       Part2(k,l)=2*sum(parameter2(k,:).*SSij_sil(k,:,l))-parameter2(k,k)*SSij_sil(k,k,l);
    end
end

Deri_of_S=Part1+Part2;


  