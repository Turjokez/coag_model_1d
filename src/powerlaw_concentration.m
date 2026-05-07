function c = powerlaw_concentration(d_cm, amp, expo)
% powerlaw_concentration
% simple power-law concentration by size

if nargin < 2 || isempty(amp)
    amp = 1.0;
end
if nargin < 3 || isempty(expo)
    expo = -2.5;
end

d_cm = d_cm(:);
ref = min(d_cm(d_cm > 0));
if isempty(ref)
    ref = 1.0;
end

c = amp .* (d_cm ./ ref) .^ expo;
c(~isfinite(c)) = 0;
c(c < 0) = 0;
end

