function [L_G, A, temp_A0, Y, distX, Wv, Parameter]=update_MUSLA_AY(XX,Parameter,Num_iter)
[mX1,nX]=size(XX{1,1});
distSUM = zeros(nX);
for ii=1:size(XX,1)
    X=XX{ii,1};
    [mX,nX]=size(X);
    distX_initial(:,:,ii) = EuDist2(X',X');
    dist = reshape(distX_initial(:,:,ii),1,nX*nX);
    Parameter.sigma1(ii) = median(dist);
    distX_initial(:,:,ii) =  (distX_initial(:,:,ii).^2)./(Parameter.sigma1(ii)^2);
    distSUM = distSUM + distX_initial(:,:,ii);
    temp_A(:,:) = exp(-distX_initial(:,:,ii));
    temp_A0{ii} = temp_A;
%     figure(ii)
%     imshow(temp_A)
end
distSUM = distSUM/size(XX,1);
k = Parameter.k_neigh;
[distX1, idx] = sort(distSUM,2);
A = zeros(nX);
for i = 1:nX
    di = distX1(i,2:k+2);
    rr(i) = 0.5*(k*di(k+1)-sum(di(1:k)));
    id = idx(i,2:k+2);
    A(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
end
r = mean(rr);
Parameter.r = r;
A = (A+A')/2;
D = diag(sum(A));
L = D - A;
[Y, temp, evs]=eig1(L, Parameter.C, 0);

Y = Y';
for iter = 1:Num_iter
    
    distSUM = zeros(nX);
    for ii=1:size(XX,1)
        if iter ==1
            distX_updated = distX_initial;
        end
        Wv(ii) = 0.5/sqrt(sum(sum( distX_updated(:,:,ii).*A)));  
        distX_updated(:,:,ii) = Wv(ii)*distX_updated(:,:,ii) ;
        distSUM = distSUM + distX_updated(:,:,ii);
    end
    distX = distSUM;
    
    if iter ==1 && Parameter.Y_is_initialization == 0 
        Y = 0*Y;
    end
    distf = EuDist2(Y',Y');
    distf = distf .^ 2;

    A = zeros(nX);
    for i=1:nX
        if Parameter.islocal == 1
            idxa0 = idx(i,2:k+1);
        else
            idxa0 = 1:num;
        end
        dfi = distf(i,idxa0);
        dxi = distX(i,idxa0);
        ad = -(dxi+Parameter.lambda_0*dfi)/(2*r);
        A(i,idxa0) = EProjSimplex_new(ad);
    end

    A = (A+A')/2;
    D_A = diag(sum(A));
    L_G=D_A-A; %Laplacian matrix of A;Spectral analysis;
    Y_old = Y ;
    [Y, temp, ev]=eig1(L_G, Parameter.C, 0);
    Y = Y';
    
    if Parameter.lambda_0_determined == 1
        fn1 = sum(ev(1:Parameter.C));
        fn2 = sum(ev(1:Parameter.C+1));
        if fn1 > 0.00000000001
            Parameter.lambda_0 = 2*Parameter.lambda_0;
        elseif fn2 < 0.00000000001
            Parameter.lambda_0 = Parameter.lambda_0/2;  Y = Y_old;
        else
            break;
        end
    end
end