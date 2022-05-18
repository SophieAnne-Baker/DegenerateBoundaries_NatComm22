dbstop if error
close all; clear all; clf
%% general to all subfigs

% > thresholds

% functions to rotate 3D probabilities to 2D plane
proj2Dpoint = @(p1,p2,p3) cat(3, sqrt(1/2)*squeeze(p3-p1+1), sqrt(2/3)*squeeze(p2-p1/2+1/2-p3/2));
proj2D = @(p) squeeze( proj2Dpoint( p(:,1)', p(:,2)', p(:,3)') );

% functions to find edge of decision boundaries
edgeDec = @(pdec, d) circshift( [pdec 1-pdec 0; pdec 0 1-pdec], [1 d-1]);
lineDec = @(p, n) p(1,:) + linspace(0,1,n)'*diff(p);

% vertices of region
vertices3D = [1 0 0; 0 1 0; 0 0 1; 1 0 0];
vertices2D = proj2D( vertices3D );

% set up meshgrid of probabilities
n = 10000;
ngrid = 10000;
[p(:,:,1), p(:,:,2)] = meshgrid( linspace( 0, 1, ngrid ) ); p(:,:,3) = 1 - sum(p,3);
p = reshape(p,[ngrid^2 3]);
p = p( p(:,3)>0 , : );

% curve boundaries 
% scaling function
scale = @(p, alpha) 1 + (1/nchoosek(size(p, 2), 3)) .* alpha .* sum( prod( reshape(p(:, nchoosek(1:size(p, 2), 3)), [size(p, 1), nchoosek(size(p, 2), 3), 3]), 3), 2);

% decision boundaries
alpha = [15.454545454545453  -8.181818181818180  -0.909090909090907]; %curve2

% find decision boundaries
pdec = [0.593591794935964   0.444818987505391   0.459876838497501]; %curve2
ndec = length(pdec);
for idec = 1 : ndec
    for d = 1:3
        boundary3D = lineDec( edgeDec( pdec(idec), d ), n );
        boundary2D{idec, d} = proj2D( boundary3D );
    end
    
    logp = log(p);
    logpdec = log( pdec(idec) );
    isdec = 0.00001 > abs( logp' - scale(p, alpha(idec))' .* ones(3,1) * logpdec );
    p2D_curve{idec} = proj2D( p( any(isdec), : ) );
end

% power boundaries
%scaling function
scale = @(p, alpha, beta) 1 + (1/nchoosek(size(p, 2), 3)) .* alpha .* sum( (max((reshape(p(:, nchoosek(1:size(p, 2), 3)), [size(p, 1), nchoosek(size(p, 2), 3), 3])), [], 3) - min((reshape(p(:, nchoosek(1:size(p, 2), 3)), [size(p, 1), nchoosek(size(p, 2), 3), 3])), [], 3)).^(beta) .* prod( reshape(p(:, nchoosek(1:size(p, 2), 3)), [size(p, 1), nchoosek(size(p, 2), 3), 3]), 3), 2);

% decision boundaries
alpha = [-18.3673469387755   -2.0408163265306 -100.0000]; %power2
beta = [0.408163265306122  -2.040816326530614   2.040816326530610]; %power2

% find decision boundaries
pdec = [0.333333333333333   0.333333333333333   0.389976937586229]; %power2
ndec = length(pdec);
for idec = 1 : ndec
    for d = 1:3
        boundary3D = lineDec( edgeDec( pdec(idec), d ), n );
        boundary2D{idec, d} = proj2D( boundary3D );
    end
    
    logp = log(p);
    logpdec = log( pdec(idec) );
    isdec = 0.00001 > abs( logp' - scale(p, alpha(idec), beta(idec))' .* ones(3,1) * logpdec );
    p2D_power{idec} = proj2D( p( any(isdec), : ) );
end

% oscil boundaries
% scaling function
scale = @(p, alpha, beta) 1 + (1/nchoosek(size(p, 2), 3)) .* alpha .* sum( cos( beta .*  2.*pi .* (max((reshape(p(:, nchoosek(1:size(p, 2), 3)), [size(p, 1), nchoosek(size(p, 2), 3), 3])), [], 3) - min((reshape(p(:, nchoosek(1:size(p, 2), 3)), [size(p, 1), nchoosek(size(p, 2), 3), 3])), [], 3)) ).* prod( reshape(p(:, nchoosek(1:size(p, 2), 3)), [size(p, 1), nchoosek(size(p, 2), 3), 3]), 3), 2); %alpha = -5.5:4.5

% decision boundaries
 alpha = [-0.612244897959183 -10.000000000000000   7.959183673469386]; %-10 %oscil3
 beta = [50.000000000000000   1.020408163265306   4.081632653061225]; %19.387755102040817 %oscil3

% find decision boundaries
pdec = [0.570914125497677   0.597097324713625   0.610635844327370]; %0.597097324713625 %oscil3
ndec = length(pdec);
for idec = 1 : ndec
    for d = 1:3
        boundary3D = lineDec( edgeDec( pdec(idec), d ), n );
        boundary2D{idec, d} = proj2D( boundary3D );
    end
    
    logp = log(p);
    logpdec = log( pdec(idec) );
    isdec = 0.00001 > abs( logp' - scale(p, alpha(idec), beta(idec))' .* ones(3,1) * logpdec );
    p2D_oscil{idec} = proj2D( p( any(isdec), : ) );
end

% > reaction times

N = 3;

% distribution
dist.generate = @generateLogp;
dist.mu = 1/2; dist.si = 1; dist.N = N;
dist.ntrials = 100000;
dist.maxt = 30;
 
% Analysis
nsamples = 500;
 
% generate distribution data
[dist.z,dist.target] = dist.generate(dist);

% symmetric -- non-target specific rTs
for i=1:nsamples
    traj2Dsym(:,:,i) = proj2D(exp(dist.z(:,:,i)));
end

clear p;
%% 

pdec = 0.35:0.01-0.005/2:1;
ndec = length(pdec);
max_logp = log(max(exp(dist.z), [], 2));
for idec = 1 : ndec
        logpdec = log( pdec(idec) );
        predec = max_logp > logpdec+0.00001;
        
        for samp = 1:size(dist.z,3)
            ind = find(predec(:,1,samp),1,'first');
            distz_bound{idec}(:,:,samp) = [dist.z(1:ind-1, :, samp); nan([32-ind+1,3, 1])];
            distz_bound_2D{idec}(:,:,samp) = proj2D(exp(distz_bound{idec}(:,:,samp)));
        end
end

cmap = colormap(parula(32));


%% 

pdec = 0.35:0.01:1;
alphdec = [];
ndec = length(pdec);
max_logp = log(max(exp(dist.z), [], 2));
for idec = 1 : ndec
        logpdec = log( pdec(idec) );
        predec = max_logp > logpdec+0.00001;
        
        for samp = 1:size(dist.z,3)
            ind = find(predec(:,1,samp),1,'first');
            distz_bound{idec}(:,:,samp) = [dist.z(1:ind-1, :, samp); nan([32-ind+1,3, 1])];
            distz_bound_2D{idec}(:,:,samp) = proj2D(exp(distz_bound{idec}(:,:,samp)));
        end
end

cmap = colormap(fliplr(turbo(32)));

cmap_temp = colormap(fliplr(turbo(96)));
cmap = ([cmap_temp(1:6:32,:); cmap_temp(33:3:64,:); cmap_temp(65:2:end-2,:)]);


%% curve 
set(0,'DefaultFigureColor',[1 1 1])
figure; hold on;
for idec = ndec:-1:1
    for t = 1:32
        scatter(distz_bound_2D{idec}(t,1,:), distz_bound_2D{idec}(t,2,:), 2, cmap(t,:));
    end
end

% plot example thresholds
hold on;
col = [0 0 1; [220,20,60]./255; 1 1 0]; % blue, red, yellow
for idec = 1:3
    plot(p2D_curve{idec}(:,1), p2D_curve{idec}(:,2), '.', 'color',col(idec,:),'markersize',7); 
end
plot(vertices2D(:,1), vertices2D(:,2), 'k')
axis off
 ha = gca; 
 ha.YLim(1)=0;
 ha.Position = ha.Position.*[1 1 0.9 1];
 set(0,'DefaultAxesTitleFontWeight','normal');
 title('curve (\theta, \alpha)', 'Units', 'normalized', 'Fontsize', 36, 'Position', [0.5, -0.12, 0]);

hc = colorbar; 
hc.Ticks = 0.1:0.2:0.9;
hc.TickLabels = 0.5:0.1:0.9; 
hc.Position = hc.Position.*[1.35 1 0.25 0.85];
hc.Label.String = 'threshold, \Theta';

if ~exist('figs', 'dir')
       mkdir('figs')
end

drawnow; pause(0.05);    

export_fig(append('./figs/', mfilename, 'a'), '-pdf', '-dpdf', '-opengl', '-eps', '-q101');
close all;

%% power
set(0,'DefaultFigureColor',[1 1 1])
figure; hold on;
for idec = ndec:-1:1
    for t = 1:32
        scatter(distz_bound_2D{idec}(t,1,:), distz_bound_2D{idec}(t,2,:), 2, cmap(t,:));
    end
end

% plot example thresholds
hold on;
col = [0 0 1; [220,20,60]./255; 1 1 0];
for idec = 1:3
    plot(p2D_power{idec}(:,1), p2D_power{idec}(:,2), '.', 'color',col(idec,:),'markersize',7); 
end
plot(vertices2D(:,1), vertices2D(:,2), 'k')
axis off
 ha = gca; 
 ha.YLim(1)=0;
 ha.Position = ha.Position.*[1 1 0.9 1];
 set(0,'DefaultAxesTitleFontWeight','normal');
 title('power (\theta, \alpha, \beta)', 'Units', 'normalized', 'Fontsize', 36, 'Position', [0.5, -0.12, 0]);

hc = colorbar; 
hc.Ticks = 0.1:0.2:0.9;
hc.TickLabels = 0.5:0.1:0.9; 
hc.Position = hc.Position.*[1.35 1 0.25 0.85];
hc.Label.String = 'threshold, \Theta';

if ~exist('figs', 'dir')
       mkdir('figs')
    end

export_fig(append('./figs/', mfilename, 'b'), '-pdf', '-dpdf', '-opengl', '-eps', '-q101');
close all;

%% oscil
set(0,'DefaultFigureColor',[1 1 1])
figure; hold on;
for idec = ndec:-1:1
    for t = 1:32
        scatter(distz_bound_2D{idec}(t,1,:), distz_bound_2D{idec}(t,2,:), 2, cmap(t,:));
    end
end

hold on;
col = [0 0 1; [220,20,60]./255; 1 1 0];
for idec = 1:3
    plot(p2D_oscil{idec}(:,1), p2D_oscil{idec}(:,2), '.', 'color',col(idec,:),'markersize',7); 
end
plot(vertices2D(:,1), vertices2D(:,2), 'k')
axis off
 ha = gca; 
 ha.YLim(1)=0;
 ax_sz = ha.Position.*[1 1 0.9 1];
 ha.Position = ax_sz;
 set(0,'DefaultAxesTitleFontWeight','normal');
 title('oscil (\theta, \alpha, \beta)', 'Units', 'normalized', 'Fontsize', 36, 'Position', [0.5, -0.12, 0]);
 
cmap = colormap(fliplr(turbo(32)));

cmap_temp = colormap(fliplr(turbo(96)));
cmap = ([cmap_temp(1:6:32,:); cmap_temp(33:3:64,:); cmap_temp(65:2:end-2,:)]);

hc = colorbar('Ticks',[0:1/6:1],'TickLabels',{'0','5','10','15','20','25','30'},'FontSize',18);
hc.Label.String = 'mean decision time';
hc.Position(4) = 0.7;

set(ha, 'Position', ax_sz);
hc.Position(4) = 0.7;

if ~exist('figs', 'dir')
       mkdir('figs')
    end

export_fig(append('./figs/', mfilename, 'c'), '-pdf', '-dpdf', '-opengl', '-eps', '-q101');
close all