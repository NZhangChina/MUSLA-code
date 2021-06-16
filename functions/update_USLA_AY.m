function [L_G A Y distX Parameter]=update_USLA_AY(X,Parameter,Num_iter)
[mX,nX]=size(X);


distX = EuDist2(X',X');
dist = reshape(distX,1,nX*nX);
Parameter.sigma1 = median(dist);
distX =  (distX.^2)./(Parameter.sigma1^2+eps);
k = Parameter.k_neigh;
[distX1, idx] = sort(distX,2);
A = zeros(nX);
for i = 1:nX
    di = distX1(i,2:k+2);
    rr(i) = 0.5*(k*di(k+1)-sum(di(1:k)));
    id = idx(i,2:k+2);
    A(i,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
end
r = mean(rr);
Parameter.r = r;
% if Parameter.lambda_0_determined == 1
%     Parameter.lambda_0 = r;
% end
A0 = (A+A')/2;
D0 = diag(sum(A0));
L0 = D0 - A0;
[Y, temp, evs]=eig1(L0, Parameter.C, 0);

Y = Y';
for iter = 1:Num_iter
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
    evs(:,iter+1) = ev;
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