function [L_G G Parameter]=Spectral_timeseries_similarity(X,Parameter)
[mX,nX]=size(X);

gg = EuDist2(X',X');
dist = reshape(gg,1,nX*nX);
Parameter.sigma1 = median(dist)+10^-10;
G = exp(- (gg.^2)./(Parameter.sigma1^2));
D_G = diag(sum(G));
L_G=D_G-G; %Laplacian matrix of G;Spectral analysis;
