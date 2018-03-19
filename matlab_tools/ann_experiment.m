function [model,results] = ann_experiment(X,Y,I,J,trnfctn,Nepochs)

% make sure no missing values
assert(isempty(find(isnan(X(:)))));
assert(isempty(find(isnan(Y(:)))));
assert(isempty(find(X(:)<-9990)));
assert(isempty(find(Y(:)<-9990)));

% train
fprintf('\nTraining: Ntrn = %d - Nepochs %d ... ',length(I),Nepochs); tic;
model = train_ann(X,Y,I,trnfctn,Nepochs);
t = toc; fprintf('finished - time = %f \n\n',t);

% test
fprintf('Predicting: Ndata = %d ... ',length(Y)); tic;
results = pred_ann(X,Y,model,I,J);
t = toc; fprintf('finished - time = %f \n\n',t);

% screen report
fprintf('Ntrn = %d    Dx = %d \n',length(I),size(X,2));
fprintf(' - Info (1%%) : %f, %f, %f \n',squeeze(results.info(:,1)));
fprintf(' - Info (2%%) : %f, %f, %f \n',squeeze(results.info(:,2)));
fprintf(' - Info (5%%) : %f, %f, %f \n',squeeze(results.info(:,3)));
fprintf(' - Info (10%%): %f, %f, %f \n',squeeze(results.info(:,4)));
fprintf(' - RMSE      : %f, %f, %f \n' ,squeeze(results.rmse));
fprintf(' - CORR      : %f, %f, %f \n' ,squeeze(results.corr));
fprintf('\n');


