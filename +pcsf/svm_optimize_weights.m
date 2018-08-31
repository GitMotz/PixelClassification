function [y] = svm_optimize_weights(x,XTraining,YTraining,TestLabel,kernel)
for i = 1:size(XTraining,2)
    XTraining(:,i)=XTraining(:,i).*x(i);
end

for i = 1:6
TestLabel = rand(length(YTraining),1)>0.33;

Model = pccore.svm_train(XTraining(TestLabel,:),YTraining(TestLabel,:),kernel);
[predictions] = predict(Model, XTraining(~TestLabel,:));
y(1,i)=1-[sum(YTraining(~TestLabel,1)==predictions)./length(predictions)];

% Model = pccore.svm_train(XTraining(~TestLabel,:),YTraining(~TestLabel,:),kernel);
% [predictions] = predict(Model, XTraining(TestLabel,:));
% y(2,i)=1-[sum(YTraining(TestLabel,1)==predictions)./length(predictions)];
end
y=mean(y(:));
