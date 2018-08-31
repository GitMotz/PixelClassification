function [Model,CModel,oberror] = svm_train(XTraining,YTraining,kernel)

%index the daata... (turned off by Nico)
%XTraining = nanzscore(XTraining,[],2);
randix = rand(size(YTraining))<0.750;%true(size(YTraining));%
%XTraining(:,17)=XTraining(:,17).*0.2;
%: = [1:size(XTraining,2)];
%XTraining = nanzscore(XTraining);
% x=[1.4325    0.4232    0.2765    1.7440    0.8761    1.8778    1.8810    1.6499    0.7314    1.8414    1.4027 ...
%    1.7162    0.7840    1.0894    1.9904    1.9953    0.4880    1.6219    0.2951    0.3065    0.3500    1.7442];
% x = [1.0676    1.9348    0.8916    0.8840    1.8055    0.6512    0.5698    1.4286    1.6081    1.8982    1.0537...
% 1.4383    0.2437    0.6411    1.8387    1.2553    0.2234    1.4112    0.2171    0.3805    1.9843    1.8130];
% for i = 1:size(XTraining,2)
%     XTraining(:,i)=XTraining(:,i).*x(i);
% end



%pre svm
% SVMStruct1 =  fitcsvm(XTraining(~randix,:),YTraining(~randix,1),...
%     'KernelFunction','Linear');
% CompactSVMModel1 = compact(SVMStruct1);
% CompactSVMModel1 = fitPosterior(CompactSVMModel1,...
%     XTraining(randix,:),YTraining(randix,1));
% [~,p] = predict(SVMStruct1, XTraining);
% XTraining = [XTraining p];

switch kernel
    case 'SVMEnsemble'
        nsvm = 1;
        [Model,CModel] = pccore.svm_ensemble_train(XTraining(randix,:),YTraining(randix,1),nsvm);
        [a,b] = pccore.svm_ensemble_predict(CModel,XTraining(~randix,:));
        predictions = mean(b,3);
        predictions = (predictions(:,1)<0.5)+1;
        oberror=1-[sum(YTraining(~randix,1)==predictions)./length(predictions)];
        
    case 'RandomForest'
        BaggedEnsemble = TreeBagger(20,XTraining(randix,:),YTraining(randix,1),...
            'Method','classification');
        [predictions,p] = predict(BaggedEnsemble, XTraining(~randix,:));
        predictions = cellfun(@str2num,predictions);
        oberror=1-[sum(YTraining(~randix,1)==predictions)./length(predictions)];
        Model = BaggedEnsemble;
        CModel = BaggedEnsemble;
        min(p(:))
        
    case 'GentleBoost'
        fprintf('%s: Calculating gentle boosting algorith, it will take a moment.\n',mfilename);
        cost = [0 10;100 0];
        Ensemble = fitensemble(XTraining(randix,:),YTraining(randix,1),'LogitBoost',400,...
            'Tree');%,'cost',cost);                %LogitBoost
        predictions = predict(Ensemble, XTraining(~randix,:));
        oberror=1-[sum(YTraining(~randix,1)==predictions)./length(predictions)];
        Model = Ensemble;
        CModel = Ensemble;
        
    otherwise
        
        good_flag = false;
        n = 0;
        
        while ~good_flag
            
            % train the svm
            fprintf('%s: the training set has %d entries.\n',mfilename,size(XTraining,1));
            fprintf('%s: Training SVM using a %s kernel.\n',mfilename,kernel);
            SVMStruct =  fitcsvm(XTraining(randix,:),YTraining(randix,1),...
                'KernelFunction',kernel,'KernelScale','auto');
            
            
            % try to compute the posterior probability. Note: this sometimes fails.
            CompactSVMModel = compact(SVMStruct);
            try
                CompactSVMModel = fitPosterior(CompactSVMModel,...
                    XTraining(randix,:),YTraining(randix,1));
            catch
                warning('Cannot fit posterior! Try another kernel.')
            end
            
            % calculate the out-of-box error
            [predictions] = predict(SVMStruct, XTraining(~randix,:));
            oberror=1-[sum(YTraining(~randix,1)==predictions)./length(predictions)];
            Model = SVMStruct;
            CModel = CompactSVMModel;
            n=n+1;
            if oberror<0.2 || n>0
                good_flag=true;
            end
            
        end
end
fprintf('%s: The out-of-box error is %.2f percent.\n',mfilename,oberror*100);

%pixelerror=YTraining(:,1)~=predictions;






% %
% figure
% for i = 1:size(XTraining,2)
%     subplot(1,3,1)
%     plot(XTraining(YTraining(:)==1,50),XTraining(YTraining(:)==1,i),'or')
%     hold on
%     plot(XTraining(YTraining(:)==2,50),XTraining(YTraining(:)==2,i),'.b')
%     plot(XTraining(pixelerror,50),XTraining(pixelerror,i),'xk')
%
%
%     hold off
%     xlim([min(XTraining(:,50)) max(XTraining(:,50))])
%     ylim([min(XTraining(:,i)) max(XTraining(:,i))])
%
%     subplot(1,3,2)
%     plot(XTraining(YTraining(:)==1,50),XTraining(YTraining(:)==1,i),'or')
%     xlim([min(XTraining(:,50)) max(XTraining(:,50))])
%     ylim([min(XTraining(:,i)) max(XTraining(:,i))])
%
%     subplot(1,3,3)
%     plot(XTraining(YTraining(:)==2,50),XTraining(YTraining(:)==2,i),'.b')
%     xlim([min(XTraining(:,50)) max(XTraining(:,50))])
%     ylim([min(XTraining(:,i)) max(XTraining(:,i))])
%
%     pause(3)
% end
%
% %
% %
%
%
% [~,PC]=princomp(XTraining);
% figure;
% plot3(PC(YTraining(:)==1,1),PC(YTraining(:)==1,2),PC(YTraining(:)==1,3),'.r')
% hold on
% plot3(PC(YTraining(:)==2,1),PC(YTraining(:)==2,2),PC(YTraining(:)==2,3),'.b')
% hold off

