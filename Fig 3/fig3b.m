
acc =  0.001;
ngrid = 1000;

tic;
% functions to rotate 3D probabilities to 2D planenotes
proj3D = @(p) p*[-1/2 -sqrt(3)/6 -sqrt(6)/12; 1/2 -sqrt(3)/6 -sqrt(6)/12; 0 sqrt(3)/3 -sqrt(6)/12; 0 0 sqrt(6)/4];

% functions to find edge of decision boundaries
    edgeDec = @(pdec, d) circshift( [pdec 1-pdec 0 0; pdec 0 1-pdec 0; pdec 0 0 1-pdec], [1 d-1]);

% vertices of region
vertices4D = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; 1 0 0 0; 1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
vertices3D = proj3D( vertices4D );

% scaling function
scale = @(p, alpha, beta) 1 + (1/nchoosek(size(p, 2), 3)) .* alpha .* sum( cos( beta .*  2.*pi .* (max((reshape(p(:, nchoosek(1:size(p, 2), 3)), [size(p, 1), nchoosek(size(p, 2), 3), 3])), [], 3) - min((reshape(p(:, nchoosek(1:size(p, 2), 3)), [size(p, 1), nchoosek(size(p, 2), 3), 3])), [], 3)) ).* prod( reshape(p(:, nchoosek(1:size(p, 2), 3)), [size(p, 1), nchoosek(size(p, 2), 3), 3]), 3), 2); %alpha = -5.5:4.5

% decision boundaries
% beta = 2;
% alpha = -10;

% beta = 2;
% alpha = 10;

beta = 5;
alpha = -5;

% set up meshgrid of probabilities
if ~exist('p', 'var')
    ngrid = 1000;
    [p(:,:,:,1), p(:,:,:,2), p(:,:,:,3)] = meshgrid( linspace( 0, 1, ngrid ) );
    p(:,:,:,4) = 1 - sum(p,4);
    p = reshape(p,[ngrid^3 4]);
    p = p( p(:,4)>0 , : );
end

% find decision boundaries
pdec = 0.5 : 0.1 : 0.9; ndec = length(pdec);
%colourings
cmap = cell(5, 3);
cmap{1, 1} = colormap(winter(64));
cmap{2, 1} = colormap(summer(64));
cmap{3, 1} = colormap(autumn(64));
cmap{4, 1} = colormap(spring(64));
cmap{5, 1} = colormap(cool(64));

for idec = 1 : ndec
    
    for d = 1:4
        boundary4D = lineDec( edgeDec( pdec(idec), d ));
        boundary3D{idec, d} = proj3D( boundary4D );
    end
    
    logp = log(p);
    logpdec = log( pdec(idec) );
    isdec = acc > abs( logp' - scale(p, alpha, beta)' .* ones(4,1) * logpdec );
    pt = p( any(isdec), : );
    ptd = pt;
    ptd(:, d) = [];
    [vals, I] = sortrows([pt, sqrt(sum(exp(ptd-1/4), 2))], [d 5], {'descend' 'ascend'});
    I = I(1:round(size(pt(:,d))/4));
    if size(unique(pt(I, d)))==1
        col{idec}.pair(1:length(I), :) = cmap{idec, 1}(32.*ones(length(I), 1), :);
    else
        col{idec}.pair(1:length(I), :) = cmap{idec, 1}((1+round( ((pt(I, d)-min(pt(I, d)))./(max(pt(I, d))-min(pt(I, d)))).*63) ), :);
    end
    for c = 0:3
        p3D{idec}(1+c*(length(I)):c*(length(I))+(length(I)), :) = proj3D( pt(I, circshift(1:1:4, c)) );
        col{idec}.pair(1+c*(length(I)):c*(length(I))+(length(I)), :) = col{idec}.pair(1:length(I), :);
    end
    
end

% plot distributions
figure(1); clf; hold on;
 colb = colormap(jet(ndec));
for idec = 1:ndec 
    scatter3(p3D{idec}(:,1), p3D{idec}(:,2), p3D{idec}(:,3), 1.5, col{idec}.pair); 
end
plot3(vertices3D(:,1), vertices3D(:,2), vertices3D(:,3), 'k')
axis off
ha = gca;

toc
set(gcf, 'color', 'w', 'units', 'inches', 'position', [0 0 2.75 2.5])
view([-17 12])


if ~exist('figs', 'dir')
       mkdir('figs')
    end

export_fig(append('./figs/', mfilename), '-pdf', '-eps', '-q101');
savefig([pwd '/figs/' mfilename])

% functions to find edge of decision boundaries
function [lD] = lineDec(p)
    step_1 = (linspace(1, 7, 7)')*reshape(diff([p;p(1,:)])', [1, 12])/8 + reshape(p', [1, 12]);
    step_2 = reshape(step_1', [4 ,7*3])';
    lD = zeros(43, 4);
    pc = circshift(p, 1, 1);
    lD(1:2:end, :) = [repmat(pc, [7, 1]); pc(end, :)];
    lD(2:2:end, :) = step_2;
    lD(end+1:end+4, :) = [p;p(1,:)];
end