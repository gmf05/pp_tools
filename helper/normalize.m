function y = normalize(x)
  mn = min(min(x)); mx = max(max(x));
  y = x-mn;
  y = y./(mx-mn);
end