function dist = kldiv(varargin)

keyboard

% extract
if isa(varargin{1},'network') || isstruct(varargin{1})
    t = varargin{2};
    y = varargin{3};
else
    t = varargin{1};
    y = varargin{2};
end

% check dimensions
t = t(:);
y = y(:);
assert(length(t) == length(y));

% % histogram bins
% Bmin = min(min(t),min(y))-1e-6;
% Bmax = max(max(t),max(y))+1e-6;
% B = Bmin:Bw:(Bmax+Bw);

M = 100;

% bin data into discrete distributions
P = histogram(t,M);
Q = histogram(y,M);
assert(length(P) == length(Q));

% normalize
Q = Q./sum(Q);
P = P./sum(P);

% calculate divergence
temp =  P.*log(P./Q);
temp(isnan(temp)) = 0; % resolving the case when P(i)==0
dist = sum(temp,2);
    

