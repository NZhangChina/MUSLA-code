function [U,nEigVec] = eig_NJW(NL1, k)

[eigVectors,eigValues] = eig(NL1);

nEigVec = eigVectors(:,(size(eigVectors,2)-(k-1)): (size(eigVectors,2))); % 样本数*标签
%nEigVec = compute_mapping(nEigVecs, 'PCA',1);

for i=1:size(nEigVec,1)
     n = sqrt(sum(nEigVec(i,:).^2));    
     U(i,:) = nEigVec(i,:) ./ n; 
end
