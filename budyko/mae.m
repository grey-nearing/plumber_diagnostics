function r = mae(X,Y,norm)

[N,D] = size(X);
assert(N == size(Y,1));
assert(D == size(Y,2));


r = mean(abs(X-Y));
return

for d = 1:D
    r(d) = mean(abs(X(:,d)-Y(:,d)));
    if     nargin == 3 && norm == 1
        r(d) = r(d) / mean(X(:,d));
    elseif nargin == 3 && norm == 2
        r(d) = r(d) / mean(Y(:,d));
    end
end


