function p = MYchitest(x1, x2)
edges = [-0.5:1:7];
a = histogram(x1, edges);
% edges = a.BinEdges;
a = a.Values;

b = histogram(x2, edges);
b = b.Values;
df = length(a) - 1;
ChiSqStatistic = sum((a-b).^2 ./ b);
p = chi2cdf(ChiSqStatistic,df,'upper')
end