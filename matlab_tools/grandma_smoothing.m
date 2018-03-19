function [X,M] = grandma_smoothing(X,cutoff,varargin)

verbose = 1;
if nargin == 3; verbose = varargin{1}; end;

% number of data points
N = length(X);

% enough valid data?
M = length(find(isnan(X)));
%if N-M<2*L; error('not enough data'); end; 

% enough good data to try gap-filling?
if M/N > cutoff 
 if verbose; fprintf('not enough good data %d, %d \n',N,M); end
 X = zeros(size(X))/0;
 M = length(X);
 return;
end

% need at least one good data to start
if isnan(X(1)); 
 if verbose; fprintf('first data point is missing \n'); end
 X = zeros(size(X))/0;
 M = length(X);
 return;
end

% count consecutive missings
seq = find(isnan(X));
count = 0;
while ~isempty(seq)
 count = count+1; 
 seq = seq(find(diff(seq)==1));
end

if count > 5
 if verbose; fprintf('too many consequtive missing values: %d  \n',count); end; 
 X = zeros(size(X))/0;
 M = length(X);
 return
end

% fill in the gaps
I = find(isnan(X),1,'first');
while ~isempty(I);
 X(I) = X(I-1);
 I = find(isnan(X),1,'first');
end 



