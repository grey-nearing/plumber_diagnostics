function Yhat = ann_train_pred(X,Y,Itrn)

% ---- set up data ---------------------------------------

assert(isempty(find(isnan(X(:)))));
assert(isempty(find(isnan(Y(:)))));

[Ny,D] = size(Y); assert(D == 1);
[Nx,D] = size(X); assert(Nx == Ny); N = Nx;

% ---- set up model --------------------------------------

% training function and number of nodes
net = feedforwardnet(ceil(size(X,2)/2)+2,'trainscg');

% data pre-processing (standard normal)
net.inputs{1}.processFcns      = {'removeconstantrows','mapstd'};
net.outputs{2}.processFcns     = {'removeconstantrows'};

% training vs. validation
net.divideParam.trainRatio     = 0.75;
net.divideParam.valRatio       = 0.25;
net.divideParam.testRatio      = 0;

% performance function
net.performFcn                 = 'mse';

% number of iterations
net.trainParam.epochs          = 200;

% convergence criteria
net.trainParam.max_fail        = max(500,round(net.trainParam.epochs/10));
net.trainParam.goal            = 1e-12;

% display flags 
net.trainParam.showWindow      = 0;
net.trainParam.showCommandLine = 0;

% scale targets between [-0.8,0.8]
Ymin = min(Y); Ymax = max(Y);
Y = 1.6*(Y-Ymin)/(Ymax-Ymin)-0.8;

% ---- do the stuff ---------------------------------------

% train the network
[net,tr] = train(net,X(Itrn,:)',Y(Itrn)');

% test
Yhat = net(X'); Yhat = Yhat';
Yhat = (Yhat+0.8)/1.6 * (Ymax-Ymin)+Ymin;









