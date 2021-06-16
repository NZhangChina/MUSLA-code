function [H XS Hij_sil Parameter]=shapelet_similarity(S,Parameter)
alpha=Parameter.alpha;
[mS,nS]=size(S);
DS=S(:,2:end);
Hij_sil=zeros(mS,mS,nS-1);
for i=1:mS
    length_s=S(i,1);
    sh_s=DS(i,1:length_s); % the i-th shapelet
    for j=i:mS
        length_l=S(j,1);
        sh_l=DS(j,1:length_l); % the j-th shapelet
        [XS(i,j) XSij_si{i,j}]=distance_longseries_shortseries(sh_l,sh_s,alpha);
        XS(j,i)= XS(i,j);% calculate distance between the i-th shapelet and the j-th shapelet 
    end
end

dist = reshape(XS,1,mS*mS);
Parameter.sigma2 = median(dist)+10^-10;

for i=1:mS
    length_s=S(i,1);
    for j=i:mS
        % calculate the similarity matrix of shapelets
        H(i,j)=exp(-XS(i,j)^2/Parameter.sigma2^2); 
        XS(j,i)=XS(i,j);
        H(j,i)=H(i,j);    
        
        % calculate the derivative of H_(ij) on S_(il)
        Hij_sil(i,j,1:length_s)=H(i,j)*(-2/Parameter.sigma2^2*XS(i,j))*XSij_si{i,j};
        Hij_sil(j,i,1:length_s)=Hij_sil(i,j,1:length_s);
    end
end
  