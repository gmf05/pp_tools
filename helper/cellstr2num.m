function out = cellstr2num(in)

N = length(in);
out = zeros(1,N);
for n = 1:N, out(n) = str2num(in{n}); end

end