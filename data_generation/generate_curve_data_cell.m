function generate_curve_data_cell(N, filename, nt, ntrials, nsamples)

%instantiating cells
dist = cell(1,size(N,2));
reward = cell(1,size(N,2));
opt = cell(1,size(N,2));

r = cell(1,size(N,2));
r_std = cell(1,size(N,2));
v = cell(1,size(N,2));
v_std = cell(1,size(N,2));

% distribution
if ~exist('N', 'var'); N=3; end
if ~exist('filename', 'var'); error('You must supply a filename.'); end
if ~exist('nt', 'var'); opt{1}.nt=100; else; opt{1}.nt=nt; end
if ~exist('ntrials', 'var'); dist{1}.ntrials=100000; else; dist{1}.ntrials=ntrials; end
if ~exist('nsamples', 'var'); reward{1}.nsamples = 10000; else; reward{1}.nsamples = nsamples; end

for k=1
    tic
    dist{k}.generate = @generateLogp;
    dist{k}.mu = 1/2; dist{k}.si = 1; dist{k}.N = N(k);
    dist{k}.maxt = 300;

% reward function
reward{k}.fun = @rewardScaled;
reward{k}.wcs = linspace(0, 0.2, 1000)';

% optimization method
opt{k}.rx = [-log(dist{k}.N) -eps; -20 20]; 

%% Analyis

% generate distribution data
[dist{k}.z,dist{k}.target] = dist{k}.generate(dist{k});

% generate samples
[x1, x2] = meshgrid( linspace(0, 1, opt{k}.nt), linspace(0,1,opt{k}.nt) );
x1 = opt{k}.rx(1,1) + x1*diff(opt{k}.rx(1,:));
x2 = opt{k}.rx(2,1) + x2*diff(opt{k}.rx(2,:));

% sample rewards
sz_x1 = size(x1,1);
sz_x2 = size(x1,2);
reward_l = reward{k};
dist_l = dist{k};
parfor i = 1:sz_x1
    for j = 1:sz_x2
        [r_l(i,j,:), r_std_l(i,j,:), v_l(i,j,:,:), v_std_l(i,j,:,:)] = reward_l.fun([x1(i,j) x2(i,j)], dist_l, reward_l.wcs, reward_l.nsamples);

    end
    disp(['The index of this iteration is: ', num2str(i)])
end
r{k} = r_l;
r_std{k} = r_std_l;
v{k} = v_l;
v_std{k} = v_std_l;

toc

end

folder = fullfile(cd, '..', 'data');
if ~exist(folder, 'dir')
    mkdir(folder)
end

save(fullfile(folder, filename),'r','r_std','v','v_std','dist','reward','opt','x1','x2', '-v7.3');

end