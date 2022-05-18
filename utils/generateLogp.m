function [z,target] = generateLogp(dist)
% generate Logp decision/evidence accumulation trajectories
% up to N=10 choice-alternatives

disp('generating logL data...')

if ~isfield(dist,'mu'); dist.mu = 1; end
if ~isfield(dist,'si'); dist.si = 1; end
if ~isfield(dist,'N'); dist.N = length(dist.mu); end
if ~isfield(dist,'prior'); dist.prior = ones(1,dist.N)/dist.N; end
if ~isfield(dist,'ntrials'); dist.ntrials = 100000; end
if ~isfield(dist,'maxt'); dist.maxt = 100; end

%matrix of simplex means
ms = zeros(10, 11);
ms(1, 1) = -1/2;
ms(1, 2) = 1/2;
ms(2, 1:2) = -sqrt(3)/6;
ms(2, 3) = sqrt(3)/3;
ms(3, 1:3) = -sqrt(6)/12;
ms(3, 4) = sqrt(6)/4;
ms(4, 1:4) = -sqrt(10)/20;
ms(4, 5) = sqrt(10)/5;
ms(5, 1:5) = -sqrt(15)/30;
ms(5, 6) = sqrt(15)/6;
ms(6, 1:6) = -sqrt(21)/42;
ms(6, 7) = sqrt(21)/7;
ms(7, 1:7) = -sqrt(28)/56;
ms(7, 8) = sqrt(28)/8;
ms(8, 1:8) = -1/12;
ms(8, 9) = 2/3;
ms(9, 1:9) = -sqrt(5)/30;
ms(9, 10) = 3*sqrt(5)/10;
ms(10, 1:10) = -sqrt(55)/110;
ms(10, 11) = sqrt(55)/11;

if length(dist.mu)==1 
    if dist.N>10; disp('N not in range'); return; end
    dist.mu = (3/sqrt(3)).*dist.mu.*ms(1:dist.N-1, 1:dist.N);
end

% randomly select targets
[ndims, ntargets] = size(dist.mu);
target = randi(ntargets,[dist.ntrials 1]);

% generate evidence
z = nan(dist.maxt+2, ntargets, dist.ntrials);
for itrial = 1:dist.ntrials
    for idim = 1:ndims
        s = normrnd(dist.mu(idim,target(itrial)),dist.si,[dist.maxt 1]);
        logl(:,:,idim) = [log(dist.prior); log( normpdf(s,dist.mu(idim,:),dist.si) )];
    end
    cumlogl = [cumsum(sum(logl,3),1); -ones(1,ntargets)/eps];
    z(:,:,itrial) = cumlogl - log( sum( exp(cumlogl), 2) );
end


