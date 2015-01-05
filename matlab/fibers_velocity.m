files = dir('input/*.velocity');

index = 0;

f = figure();
set(f,'Units','Inches');
set(f, 'Position', [0 0 12 4]);
timesteps = size(files,1);
Y = []

for file = files'
    input=['input/' file.name]; 
    data = dlmread(input,'',1, 0);
    
    Y = [Y data(1,3)];
    
    clf;
    
    plot(Y);
    grid off;
    
    output=['output/velocity_' sprintf('%05d',index) '.pdf'];
    
    axis([0 200 -4.5 -1.5]);
    set(gca,'FontName','Fira Mono');
    set(gca,'YTick',-4.5:.5:1.5);
    tightInset = get(gca, 'TightInset');
    position(1) = tightInset(1);
    position(2) = tightInset(2);
    position(3) = 1 - tightInset(1) - tightInset(3);
    position(4) = 1 - tightInset(2) - tightInset(4);
    set(gca, 'Position', position);
    pos = get(f,'Position');
    set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    
    saveas(f, output);
    
    index = index + 1
    
    %close all;
end

