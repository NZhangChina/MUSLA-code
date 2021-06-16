function Fun=function_value_MUSLA_Simplify(XX,distX,Y,L_A,A,HH,Parameter)
part0=0.5*trace(distX*A);
part1=Parameter.lambda_0*trace(Y*L_A*Y');
part2=0;
for ii=1:size(XX,1)
    H=HH{ii,1};
    part2=part2+0.5*Parameter.lambda_1(ii)*trace(H'*H);
end
Fun=part0+part1+part2;