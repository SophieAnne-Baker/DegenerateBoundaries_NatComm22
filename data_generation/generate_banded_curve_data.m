function generate_banded_curve_data(N, filename, ntrials, nsamples)

% distribution
if ~exist('N', 'var'); N=3; end
if ~exist('filename', 'var'); error('You must supply a filename.'); end
if ~exist('nt', 'var'); opt.nt=2500; else; opt.nt=nt; end
if ~exist('ntrials', 'var'); dist.ntrials=100000; else; dist.ntrials=ntrials; end
if ~exist('nsamples', 'var'); reward.nsamples = 10000; else; reward.nsamples = nsamples; end

for k=1
    tic
    dist.generate = @generateLogp; 
    dist.mu = 1/2; dist.si = 1; dist.N = N(k);
    dist.maxt = 300;

% reward function
reward.fun = @rewardScaled;
reward.wcs = linspace(0, 0.2, 1000)';

% optimization method
opt.rx = [-log(dist.N) -eps; -20 20]; 

%% Analyis

% generate distribution data
[dist.z,dist.target] = dist.generate(dist);

% generate samples
[x1, x2] = meshgrid( linspace(0, 1, opt.nt), linspace(0,1,5) );
x1 = opt.rx(1,1) + x1*diff(opt.rx(1,:));
x2 = opt.rx(2,1) + x2*diff(opt.rx(2,:));

% sample rewards
sz_x1 = size(x1,1);
sz_x2 = size(x1,2);
reward_l = reward;
dist_l = dist;
parfor i = 1:sz_x1
    for j = 1:sz_x2
        [r_l(i,j,:), r_std_l(i,j,:), v_l(i,j,:,:), v_std_l(i,j,:,:)] = reward_l.fun([x1(i,j) x2(i,j)], dist_l, reward_l.wcs, reward_l.nsamples);

    end
    disp(['The index of this iteration is: ', num2str(i)])
end
r = r_l;
r_std = r_std_l;
v = v_l;
v_std = v_std_l;

toc

end

folder = fullfile(cd, '..', 'data');
if ~exist(folder, 'dir')
    mkdir(folder)
end

save(fullfile(folder, filename),'r','r_std','v','v_std','dist','reward','opt','x1','x2', '-v7.3');

end