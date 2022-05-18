function [r, r_std, v, v_std] = rewardScaled(x, dist, wcs, nsamples, itrial)
% x = thresholds
% wcs = cost parameters
% nsamples = # samples to average reward over

if ~exist('wcs','var') || isempty(wcs); wcs = 0.02; end
if ~exist('nsamples','var') || isempty(nsamples); nsamples = 1; end
if ~exist('itrial','var') || isempty(itrial); itrial = randi(dist.ntrials,[nsamples 1]); end

% if only one input
if isfield(dist,'N'); N = dist.N; else N = size(dist.z,3); end 
if size(wcs, 2)==1; wcs = wcs*ones(1,N); end
if size(wcs, 2)~=N; disp('2nd dimension of wcs must equal either 1 or N'); end
if length(x)==1; x(2) = eps; end

% scaling
scale = @(p, alpha) 1 + (1/nchoosek(size(p, 2), 3)) .* alpha .* sum( prod( reshape(p(:, nchoosek(1:size(p, 2), 3)), [size(p, 1), nchoosek(size(p, 2), 3), 3]), 3), 2);

% normalize weights
c = max(wcs, [], 2); ws = wcs./c;

% threshold crossing
[e, rt] = deal(nan(nsamples,1));
for isample = 1:nsamples
    sc = scale( exp( dist.z( 2:end, :, itrial(isample) ) ), x(2) );
    sc = [sc(1:end-1); 0]; % solve crossing problem
    isdec = dist.z(2:end, :, itrial(isample)) > x(1)*sc.*ones(1,N);
    rt(isample) = find( sum( isdec ,2), 1);
    [~, n] = max( dist.z(rt(isample)+1, :, itrial(isample)) );
    e(isample) = n~=dist.target(itrial(isample));
end
rnan = rt>dist.maxt; e(rnan) = nan; rt(rnan) = nan;
r = nanmean( ( - c*rt' - ws(:,dist.target(itrial)).*repmat(e', [size(wcs, 1), 1]) ), 2);
r_std = nanstd( ( - c*rt' - ws(:,dist.target(itrial)).*repmat(e', [size(wcs, 1), 1]) ),[], 2);
v = [nanmean(e) nanmean(rt)]; 
v_std = [nanstd(e) nanstd(rt)];