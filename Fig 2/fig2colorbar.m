clear all; clf;
col = colormap(jet(5));
figure(1);
set(gca,'visible','off')
c = colorbar('FontSize',16,'Ticks', [0.1, 0.3, 0.5, 0.7, 0.9], 'TickLabels',{'0.5','0.6','0.7','0.8','0.9'});
c.Label.String = 'threshold, \theta';
c.Position = [0.5    0.1096    0.04    0.8148]; %[x1 y1 width height]
set(gcf, 'color', 'w', 'units', 'inches', 'position', [0 0 2 4.5])

if ~exist('figs', 'dir')
       mkdir('figs')
    end

export_fig(append('./figs/', mfilename), '-pdf', '-eps', '-q101');
savefig([pwd '/figs/' mfilename])