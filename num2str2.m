function s = num2str2(x, p, el, eu)
% s = num2str2(x, p, el, eu). Converts x to string in exponential form when 
% below el, p precision sigfigs from el to 10, rounded integer from 10 to 
% eu, and exponential above eu
%
% Mandatory input
%   x -> Number to be rounded


if ~exist('p', 'var') || isempty(p)
    p = 2;
end
if ~exist('el', 'var') || isempty(el)
    el = 0.01;
end
if ~exist('eu', 'var') || isempty(eu)
    eu = 10000;
end

x_abs = abs(x);
if x_abs < el || x_abs > eu
    s = sprintf('%.*e', p, x);
elseif x_abs >= 10
    s = sprintf('%.0f', x);
else
    s = sprintf('%.*g', p, x);
end