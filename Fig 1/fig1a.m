clear all; clf; %close all
dbstop if error
rng(1);

nTrials = 2000;
tMax = 50;

% decision threshold
pdec = 0.66;

% parameters for distribution
dist.mu = 1/2; dist.si = 1; dist.N = 3;
dist.ntrials = nTrials;
dist.maxt = tMax;
dist.generate = @generateLogp;

%% analysis

% generate distribution
[dist.z, dist.target] = dist.generate(dist);

% decisions via threshold crossing
r = nan(1,nTrials); 
v = nan(2,nTrials);
for itrial = 1:nTrials
    [r(itrial), v(:,itrial)] = rewardWald(log(pdec), dist, [],[], itrial); % log L
end
e = v(1,:); rt = v(2,:);
nfinish = sum(isfinite(e));
mrt = nanmean(rt);

% no trajectory after decision
for itrial = 1:nTrials
    dist.z( v(2,itrial) + 2:end, :, itrial) = nan; 
end

%% plot results

rTrials = (1:20) + 40;
if dist.N == 4
rTrials = [4, 200, 100, 420, 720];
else
    rTrials = [(1:20) + 43;];
    rTrials = [4, 200, 100, 420, 720]+54;
end

% transformation to rotate coordinates
% functions to rotate 3D probabilities to 2D plane
proj3D = @(p) p*[-1/2 -sqrt(3)/6 -sqrt(6)/12; 1/2 -sqrt(3)/6 -sqrt(6)/12; 0 sqrt(3)/3 -sqrt(6)/12; 0 0 sqrt(6)/4];


% vertices of region
vertices4D = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; 1 0 0 0; 1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
vertices3D = proj3D( vertices4D );

% decision boundaries
% functions to find edge of decision boundaries
    edgeDec = @(pdec, d) circshift( [pdec 1-pdec 0 0; pdec 0 1-pdec 0; pdec 0 0 1-pdec], [1 d-1]);
    
for d = 1:4
    boundary4D = lineDec( edgeDec( pdec, d ));
    boundary3D{d} = proj3D( boundary4D );
end


% trajectories
if dist.N == 4
    for i = 1:size(dist.z, 3)
        prot(:,:,i) = proj3D(exp(dist.z(:,:,i)));
    end
else
    for i = 1:size(dist.z, 3)
        prot(:,:,i) = proj3D( [exp(dist.z(:,:,i)), ones(size(exp(dist.z(:,:,i)), 1), 1)] );
    end
end

%% 3AFC vector view

% plot distributions
figure(1); clf

hsp = subplot(131); hsp.Position = hsp.Position.*[1 1 1.1 1.1] - [0.02 0.02 0 0]; hold on
hsp.Position = hsp.Position.*0.9;
line([1 0],[0 1],[0 0],'color','k')
line([1 0],[0 0],[0 1],'color','k')
line([0 0],[1 0],[0 1],'color','k')
patch([pdec(1) pdec(1) 1],[1-pdec(1) 0 0],[0 1-pdec(1) 0],0.9*[1 1 1],'EdgeColor','none')
patch([1-pdec(1) 0 0],[pdec(1) pdec(1) 1],[0 1-pdec(1) 0],0.9*[1 1 1],'EdgeColor','none')
patch([1-pdec(1) 0 0],[0 1-pdec(1) 0],[pdec(1) pdec(1) 1],0.9*[1 1 1],'EdgeColor','none')

C = [100, 143, 255; 120, 94, 240; 220, 38, 127; 254, 97, 0; 255, 176, 0]./255;

h = plot3(squeeze(exp(dist.z(:,1,rTrials))),squeeze(exp(dist.z(:,2,rTrials))),squeeze(exp(dist.z(:,3,rTrials))),'.-');
set(h, {'color'}, num2cell(C,2));

plot3([0 0 0 1],[1 1 0 0],[0 1 1 1],':k')
plot3([pdec(1) pdec(1)],[1 1/3],[0 0],':k')
plot3([1 1/3],[pdec(1) pdec(1)],[0 0],':k')
plot3([pdec(1) pdec(1)],[0 0],[1 1/3],':k')
plot3([1 1/3],[0 0],[pdec(1) pdec(1)],':k')
plot3([0 0],[pdec(1) pdec(1)],[1 1/3],':k')
plot3([0 0],[1 1/3],[pdec(1) pdec(1)],':k')
plot3([pdec(1) pdec(1)],[1-pdec(1) 0],[0 1-pdec(1)],'--k')
plot3([1-pdec(1) 0],[pdec(1) pdec(1)],[0 1-pdec(1)],'--k')
plot3([1-pdec(1) 0],[0 1-pdec(1)],[pdec(1) pdec(1)],'--k')
title('{\bf \fontsize{10}3-Choice}');
set(gca,'fontsize',7,'ztick',0:0.33:1,'zticklabel',{'0','1/3','\Theta','1'},...
    'xtick',0:0.33:1,'xticklabel',{'0','1/3','\Theta','1'},...
    'ytick',0:0.33:1,'yticklabel',{'0','1/3','\Theta','1'})
set(gca,'fontsize',7)
text(0.9,0.2,0.2,'choose H_0','horizontalalign','center','verticalalign','middle','fontsize',8,'rotation',-60)
text(0.1,0.1,0.9,'choose H_2','horizontalalign','center','verticalalign','middle','fontsize',8)
text(0.2,0.9,0.2,'choose H_1','horizontalalign','center','verticalalign','middle','fontsize',8,'rotation',60)

view([1,1,1])
hzl = zlabel({'{\bf \fontsize{10}Diffusion in nD}','\fontsize{8}P(H_2|x)'}); hzl.Position = hzl.Position + [-0.05 0 0];
hxl = xlabel('P(H_0|x)','fontsize',8); hxl.Position = hxl.Position + [-0.15 -0.15 0];
hyl = ylabel('P(H_1|x)','fontsize',8); hyl.Position = hyl.Position + [-0.15 -0.15 0];
plot3(1/3,1/3,1/3,'.k','markersize',5)
hAxis.YRuler.FirstCrossoverValue  = 1;
set(gcf, 'color', 'w', 'units', 'inches', 'position', [0 0 7 2])

if ~exist('figs', 'dir')
       mkdir('figs')
    end

export_fig(append('./figs/', mfilename), '-pdf', '-eps', '-q101');
savefig([pwd '/figs/' mfilename])

function [lD] = lineDec(p)
        step_1 = (linspace(1, 7, 7)')*reshape(diff([p;p(1,:)])', [1, 12])/8 + reshape(p', [1, 12]);
        step_2 = reshape(step_1', [4 ,7*3])';
        lD = zeros(43, 4);
        pc = circshift(p, 1, 1);
        lD(1:2:end, :) = [repmat(pc, [7, 1]); pc(end, :)];
        lD(2:2:end, :) = step_2;
        lD(end+1:end+4, :) = [p;p(1,:)];
end
    