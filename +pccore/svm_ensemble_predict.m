function [a,b] = svm_ensemble_predict(CModels,X)
nsvm = length(CModels);
for i = 1:nsvm

   [a(:,:,i), b(:,:,i)]  = predict(CModels(i).model, X(:,CModels(i).randfeat)); 
end




