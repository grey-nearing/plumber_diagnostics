function [Xnorm,bias,scale] = standardize(X,varargin)

% bias = 0;
% scale = 1;
% Xnorm = X;
% 
% return
% -------------------------------------------------------------------------
[N,D,M] = size(X);
% -------------------------------------------------------------------------
if M > 1
    XX = zeros(N*M,D);
    for m = 1:M
        XX(N*(m-1)+1:N*m,:) = X(:,:,m);
    end
    X = XX;
end
% -------------------------------------------------------------------------
if nargin == 1
    bias  = mean(X);                      
    scale = std(X);
    scale(scale==0) = 1;
% -------------------------------------------------------------------------
elseif nargin == 3
    bias = 0;
    scale = 0;
    if ~isempty(varargin{1})
        bias = varargin{1};
    end
    if ~isempty(varargin{2})
        scale = varargin{2};
    end
% -------------------------------------------------------------------------
else
    error('Must supply either 1 or 3 inputs');
end
% -------------------------------------------------------------------------
Xnorm = (X-repmat(bias,N*M,1))./repmat(scale,N*M,1);
% -------------------------------------------------------------------------
if M > 1
    XX = zeros(N,D,M);
    for m = 1:M
        XX(:,:,m) = Xnorm(N*(m-1)+1:N*m,:);
    end
    Xnorm = XX;
end
% -------------------------------------------------------------------------
