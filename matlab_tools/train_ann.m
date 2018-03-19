function model = train_ann(X,Y,I,trnfctn,Nepochs)

% training function and number of nodes
net = feedforwardnet(ceil(size(X,2)/2)+2,trnfctn);

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
net.trainParam.epochs          = Nepochs;

% convergence criteria
net.trainParam.max_fail        = round(500);
net.trainParam.goal            = 1e-12;

% display flags 
net.trainParam.showWindow      = 0;
net.trainParam.showCommandLine = 1;

% scale targets between [-0.8,0.8]
Ymin = min(Y); Ymax = max(Y);
Y = 1.6*(Y-Ymin)/(Ymax-Ymin)-0.8;

% train the network
[net,tr] = train(net,X(I,:)',Y(I,:)');

% unscale targets
Y = (Y+0.8)/1.6 * (Ymax-Ymin)+Ymin;

% create model structure
model.ann = net;
model.Ymin = Ymin;
model.Ymax = Ymax;

