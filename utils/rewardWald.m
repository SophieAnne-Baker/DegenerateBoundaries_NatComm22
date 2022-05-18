function [r, v, r_std, v_std] = rewardWald(x, dist, wcs, nsamples, itrial)
% x = thresholds
% wcs = cost parameters
% nsamples = # samples to average reward over

if ~exist('wcs','var') || isempty(wcs); wcs = 0.02; end
if ~exist('nsamples','var') || isempty(nsamples); nsamples = 1; end
if ~exist('itrial','var'); itrial = randi(dist.ntrials,[nsamples 1]); end

% if only one input
if isfield(dist,'N'); N = dist.N; else N = size(dist.z,3); end 
if length(x)<N; x(end+1:N) = x(end); end
if length(wcs)==1; wcs = wcs*ones(1,N); end

% normalize weights
c = max(wcs); ws = wcs/c;

% threshold crossing
[e, rt] = deal(nan(nsamples,1)); 
for isample = 1:nsamples
    rt(isample) = find( sum( dist.z(2:end, :, itrial(isample))' > x(:) ,1), 1);
    [~, n] = max( dist.z(rt(isample)+1, :, itrial(isample)) ); 
    e(isample) = n~=dist.target(itrial(isample));
end
rnan = rt>dist.maxt; e(rnan) = nan; rt(rnan) = nan;

% output means
r = nanmean( - c*rt - ws(dist.target(itrial))'.*e );
r_std = nanstd( - c*rt - ws(dist.target(itrial))'.*e );
v = [nanmean(e) nanmean(rt)]; 
v_std = [nanstd(e) nanstd(rt)]; 
