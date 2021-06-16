function W_tp1=update_W(X_t,Y_t,Parameter)
[mX,nX]=size(X_t);
P1=Parameter.lambda_2*X_t*X_t'+Parameter.lambda_3*eye(mX);
P2=Parameter.lambda_2* X_t* Y_t';
W_tp1=inv(P1)*P2;