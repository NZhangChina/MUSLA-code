function Fun=function_value_USLA_Simplify(distX,Y,L_A,A,H,Parameter)
part0=0.5*trace(distX*A);
part1=Parameter.lambda_0*trace(Y*L_A*Y');
part2=0.5*Parameter.lambda_1*trace(H'*H);
Fun=part0+part1+part2;