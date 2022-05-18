clear all; clf; close all

tic;
% functions to rotate 3D probabilities to 2D plane
proj2Dpoint = @(p1,p2,p3) cat(3, sqrt(1/2)*squeeze(p3-p1+1), sqrt(2/3)*squeeze(p2-p1/2+1/2-p3/2));
proj2D = @(p) squeeze( proj2Dpoint( p(:,1)', p(:,2)', p(:,3)') );

% functions to find edge of decision boundaries
edgeDec = @(pdec, d) circshift( [pdec 1-pdec 0; pdec 0 1-pdec], [1 d-1]);
lineDec = @(p, n) p(1,:) + linspace(0,1,n)'*diff(p);

% vertices of region
vertices3D = [1 0 0; 0 1 0; 0 0 1; 1 0 0];
vertices2D = proj2D( vertices3D );

% scaling function
scale = @(p, alpha) 1 + (1/nchoosek(size(p, 2), 3)) .* alpha .* sum( prod( reshape(p(:, nchoosek(1:size(p, 2), 3)), [size(p, 1), nchoosek(size(p, 2), 3), 3]), 3), 2);

% decision boundaries
n = 10000;
alpha = 10;

% set up meshgrid of probabilities
ngrid = 10000;
[p(:,:,1), p(:,:,2)] = meshgrid( linspace( 0, 1, ngrid ) ); p(:,:,3) = 1 - sum(p,3);
p = reshape(p,[ngrid^2 3]);
p = p( p(:,3)>0 , : );

% find decision boundaries
pdec = 0.5 : 0.1 : 0.9; ndec = length(pdec);
for idec = 1 : ndec
    for d = 1:3
        boundary3D = lineDec( edgeDec( pdec(idec), d ), n );
        boundary2D{idec, d} = proj2D( boundary3D );
    end
    
    logp = log(p);
    logpdec = log( pdec(idec) );
    isdec = 0.00001 > abs( logp' - scale(p, alpha)' .* ones(3,1) * logpdec );
    p2D{idec} = proj2D( p( any(isdec), : ) );
end

% plot distributions
figure(1); clf; hold on;
col = colormap(jet(ndec));
for idec = 1:ndec
    plot(p2D{idec}(:,1), p2D{idec}(:,2), '.', 'color',col(idec,:),'markersize',1); 
    for d = 1:3 
        plot(boundary2D{idec,d}(:,1), boundary2D{idec,d}(:,2),...
            '--', 'color',col(idec,:), 'linewidth', 1); 
    end
end
plot(vertices2D(:,1), vertices2D(:,2), 'k')
axis off
text(vertices2D(1,1), 0.66*vertices2D(2,2), ['\alpha=' int2str(alpha)],'fontsize',10)

toc
set(gcf, 'color', 'w', 'units', 'inches', 'position', [0 0 2.75 2.5])

if ~exist('figs', 'dir')
       mkdir('figs')
    end

export_fig(append('./figs/', mfilename), '-pdf', '-eps', '-q101');
savefig([pwd '/figs/' mfilename])