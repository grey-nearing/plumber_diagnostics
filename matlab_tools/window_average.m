function Xw = window_average(X,W)

% get dimensions
[N,D] = size(X);
assert(N > W*2);

% series length
L = floor(N/W)-1;

% init storage
Xw = zeros(L,W)./0;

% loop through dimensions of input data
for w = 1:W
 edex = w-1;
 for l = 1:L
  sdex = edex+1;
  edex = sdex+W-1;
%  if length(find(isnan(X(sdex:edex))))<0.1*W
   Xw(l,w) = mean(X(sdex:edex));
%  else
%   Xw = zeros(L,W)./0;
%   return
%  end
 end
end

% remove most fo the redundancy
%if W > 2
% keep = [1,round(W/2)];
% Xw = Xw(:,keep);
%end

% make sure all data is filled
assert(isempty(find(isnan(Xw),1,'first')))

