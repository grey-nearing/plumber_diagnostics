function [stats,Yhat] = pred_ann(X,Y,model,I,J)

% total vector (to deal with missings)
K = [I(:);J(:)];

% make predictions 
Yhat = model.ann(X')';

% rescale
Yhat = (Yhat+0.8)/1.6 * (model.Ymax-model.Ymin)+model.Ymin;

% bounds
Ymin = min(min(Yhat),min(Y));
Ymax = max(max(Yhat),max(Y));

% infomation bins
Bw = (model.Ymax-model.Ymin)/10;  B10 = (Ymin-Bw):Bw:(Ymax+Bw); B10 = B10(:);
Bw = (model.Ymax-model.Ymin)/20;  B5  = (Ymin-Bw):Bw:(Ymax+Bw); B5  = B5(:);
Bw = (model.Ymax-model.Ymin)/50;  B2  = (Ymin-Bw):Bw:(Ymax+Bw); B2  = B2(:);
Bw = (model.Ymax-model.Ymin)/100; B1  = (Ymin-Bw):Bw:(Ymax+Bw); B1  = B1(:);

% calculate performance metrics
RMSE(1) = rmse(Y(K),Yhat(K));
RMSE(2) = rmse(Y(I),Yhat(I));
RMSE(3) = 0/0;

CORR(1) = corr(Y(K),Yhat(K));
CORR(2) = corr(Y(I),Yhat(I));
CORR(3) = 0/0;

INFO(1,1) = info(Y(K),Yhat(K),B1,B1,1);
INFO(2,1) = info(Y(I),Yhat(I),B1,B1,1);
INFO(3,1) = 0/0;

INFO(1,2) = info(Y(K),Yhat(K),B2,B2,1);
INFO(2,2) = info(Y(I),Yhat(I),B2,B2,1);
INFO(3,2) = 0/0;

INFO(1,3) = info(Y(K),Yhat(K),B5,B5,1);
INFO(2,3) = info(Y(I),Yhat(I),B5,B5,1);
INFO(3,3) = 0/0;

INFO(1,4) = info(Y(K),Yhat(K),B10,B10,1);
INFO(2,4) = info(Y(I),Yhat(I),B10,B10,1);
INFO(3,4) = 0/0;

if length(J)>1
 RMSE(3) = rmse(Y(J),Yhat(J));
 CORR(3) = corr(Y(J),Yhat(J));
 INFO(3,1) = info(Y(J),Yhat(J),B1,B1,1);
 INFO(3,2) = info(Y(J),Yhat(J),B2,B2,1);
 INFO(3,3) = info(Y(J),Yhat(J),B5,B5,1);
 INFO(3,4) = info(Y(J),Yhat(J),B10,B10,1);
end

% store stats
stats.info = INFO;
stats.rmse = RMSE;
stats.corr = CORR;
stats.Yhat = Yhat;
stats.Yobs = Y;

