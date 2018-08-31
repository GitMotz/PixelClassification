function [Models,CModels] = svm_ensemble_train(X,Y,nsvm)
sample_number = size(X,1);
feature_number = size(X,2);

randix = ( rand(sample_number,nsvm)  > 0 ) + 1;  %randi(2,sample_number,nsvm);

for i = 1:nsvm
    randfeat(:,i) = ( rand(feature_number,1)  < (rand(1)+0.2) ) + 1 ;
end
randfeat(:,1) = 2;
CModels = struct();
n = 1;
for i = 1:nsvm
   disp(i) 
   Models(i).model =  fitcsvm(X(randix(:,i)>1,randfeat(:,i)>1),Y(randix(:,i)>1,1),...
       'KernelFunction','Gaussian','CrossVal','on','Kfold',4);
   Models(i).randfeat = randfeat(:,i)>1;
   for j=1:length(Models(i).model.Trained)
    CompactSVMModel = Models(i).model.Trained{j};%compact(Models(i).model);
    CModels(n).model = fitPosterior(CompactSVMModel,X(randix(:,i)>1,randfeat(:,i)>1),Y(randix(:,i)>1,1));
    CModels(n).randfeat = randfeat(:,i)>1;
    n=n+1;
   end
end

   
   



