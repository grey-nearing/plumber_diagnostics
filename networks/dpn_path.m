function [T,H,S] = dpn_path(X,Y,lag,Bx,By,annFlag)

% get dimensions
[Nx,D] = size(Y); assert(D==1);
[Ny,D] = size(X); assert(D==1);
assert(Ny == Nx); N = Nx;

% deal with lags
Xw = window_average(X,lag); 
Yw = window_average(Y,lag);
assert(numel(Xw) == length(Yw));
if numel(Xw)<1000; return; end;

% time shift
Xt = Xw(1:end-1,:); Xt = Xt(:);
Yt = Yw(1:end-1,:); Yt = Yt(:);
Ys = Yw(2:end,:);   Ys = Ys(:);
N = length(Ys);

% make sure we have dealt with missing values
assert(isempty(find(isnan(Xt),1,'first')));
assert(isempty(find(isnan(Yt),1,'first')));
assert(isempty(find(isnan(Ys),1,'first')));

% special case when targeting self
if max(abs(X-Y)==0)
 Yt = rand(size(Yt)) * (max(Yt)-min(Yt)) + min(Yt);
end

% indexes to calculate the stats
Itrn = 1:2:N;
Itst = 2:2:N;

% if this is an ann run
if annFlag
 Ys = ann_train_pred([Xt,Yt],Ys,Itrn);
end

% do the actaul calculations
[T,H,S] = transfer_entropy(Ys(Itst),Xt(Itst),Yt(Itst),Bx,By);


