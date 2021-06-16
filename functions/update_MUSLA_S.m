function SS_tp1=update_MUSLA_S(XX,ST_t,ST0_t,SS_t,Xkj_skl,SSij_sil,SS,Wv,Parameter)

for ii=1:size(XX,1)
    X=XX{ii,1};
    S_t=SS_t{ii,1};
    DS_t=S_t(:,2:end);
    deri_of_SS{ii,1}=derivation_of_MUSLA_S(X,ST_t,ST0_t,S_t,Xkj_skl{ii,1},SSij_sil{ii,1},SS{ii,1},Wv,Parameter,ii);
end

for ii=1:size(XX,1)
    S_t=SS_t{ii,1};
    DS_t=S_t(:,2:end);
    deri_of_S=deri_of_SS{ii,1};
    i=1;
    while i<Parameter.Imax+1 
        DS_t=S_t(:,2:end);
        DS_t=DS_t-Parameter.eta*deri_of_S;
        S_t=[S_t(:,1) DS_t];
        i=i+1;
    end
    S_tp1=S_t;
    S_tp1=[S_tp1(:,1), z_regularization(S_tp1(:,2:end))]; 
    SS_tp1{ii,1} = S_tp1;
end