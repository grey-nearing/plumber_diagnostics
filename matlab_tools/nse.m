function r = nse(O,M);

[N,D] = size(M);
assert(N == size(O,1));
assert(D == size(O,2));

for d = 1:D
 num = (M-O).^2;
 den = (O-mean(O)).^2;
 r(d) = 1-sum(num)/sum(den);
end


