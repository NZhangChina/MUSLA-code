function Y_tp1=update_Y(W_tp1,X_t,L_Gt,Parameter)
P1=Parameter.lambda_2*W_tp1'*X_t;
[mL,nL]=size(L_Gt);
P2=Parameter.lambda_0*L_Gt+Parameter.lambda_2*eye(mL);
Y_tp1_hat=P1*inv(P2);

Y_tp1=Y_tp1_hat;
