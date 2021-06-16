function Fun=function_value(X,Y,L_G,H,W,Parameter)
part1=0.5*Parameter.lambda_0*trace(Y*L_G*Y');
part2=0.5*Parameter.lambda_1*trace(H'*H);
part3=0.5*Parameter.lambda_2*trace((W'*X-Y)'*(W'*X-Y));
part4=0.5*Parameter.lambda_3*trace(W'*W);
Fun=part1+part2+part3+part4;