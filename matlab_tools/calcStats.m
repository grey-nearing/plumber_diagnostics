function stats = calcStats(O,M,Bw)

% reshape
O = O(:);
M = M(:);
assert(length(O) == length(M))

% remove gradmas
Nbefore = length(O);
Im = find(isnan(O)); O(Im) = []; M(Im) = [];
Im = find(isnan(M)); O(Im) = []; M(Im) = [];

if length(O) > 0.5*Nbefore
    
    % root mean squared error
    stats.rmse = sqrt(nanmean((M-O).^2));
    
    % mean biased error
    stats.mbe = sqrt(nanmean(abs(M-O)));
    
    % correlation coeficient
    stats.r = corr(O,M);
    
    % normalized mean error
    stats.nme = sum(abs(M-O))/sum(abs(O-mean(O)));
    
    % fifth and ninety-fifth percentiles
    stats.p5  = abs(prctile(M,5 )-prctile(O,5 ));
    stats.p95 = abs(prctile(M,95)-prctile(O,95));

    % standard deviation
    stats.sd    = abs(1-std(M)/std(O));
    
    % skewness
    stats.skew  = abs(skewness(M)-skewness(O));
    
    % kurtosis
    stats.kurt  = abs(kurtosis(M)-kurtosis(O));
    
    % overlap
    hm = histcounts(M,100)./length(M);
    ho = histcounts(O,100)./length(O);
    stats.over  = sum(min(hm,ho));
    
    % mutual information
    stats.mi = mutual_info_ratio(O,M,Bw);
    
else
    
    stats.rmse  = 0/0;
    stats.mbe   = 0/0;
    stats.r     = 0/0;
    stats.nme   = 0/0;
    stats.p5    = 0/0;
    stats.p95   = 0/0;
    stats.sd    = 0/0;
    stats.skew  = 0/0;
    stats.kurt  = 0/0;
    stats.over  = 0/0;
    stats.mi    = 0/0;
    
end
