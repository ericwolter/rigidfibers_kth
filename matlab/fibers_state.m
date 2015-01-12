% files = dir('input/*.state');
files = dir('input_wip/*.out');
% files = files(1455:end);
index = 0;
f = figure();


%r = 600; % pixels per inch
%set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 4000 25000]/r);
% % set(f,'Units','Inches');
% % set(f, 'Position', [0 0 10 17]);
% tightInset = get(gca, 'TightInset');
% position(1) = tightInset(1);
% position(2) = tightInset(2);
% position(3) = 1 - tightInset(1) - tightInset(3);
% position(4) = 1 - tightInset(2) - tightInset(4);
% set(gca, 'Position', position);
% set(f,'Units','Inches');
% pos = get(f,'Position');
% set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

for file = files'
    input=['input_wip/' file.name]; 
    data = dlmread(input,'',1, 0);

    data = reshape(data',3,2,[]);
    pos = squeeze(data(:,1,:))';
    orient = squeeze(data(:,2,:))';

    Xc = pos(:,1);
    Yc = pos(:,2);
    Zc = pos(:,3);

    Tx = orient(:,1);
    Ty = orient(:,2);
    Tz = orient(:,3);   
    Xp = [Xc'+Tx'*-1;Xc'+Tx'*0;Xc'+Tx'*1];
    Yp = [Yc'+Ty'*-1;Yc'+Ty'*0;Yc'+Ty'*1];
    Zp = [Zc'+Tz'*-1;Zc'+Tz'*0;Zc'+Tz'*1];
    
    clf;
    
    plot3(Xp,Yp,Zp,'LineWidth',0.5);
    az = 0;
    el = 90;
    %view(az, el);
    %camproj('perspective');
    grid on;
    
    output=['output_wip/state_' sprintf('%05d',index) '.pdf'];
    
    % Options for Ring Simulation
%     axis([-2 2 -2 2 -2 2]);
%     set(gca,'FontName','Fira Mono');
%     set(gca,'XTick',-2:1:2);
%     set(gca,'YTick',-2:1:2);
%     zm = Zp(2,1);
%     zlim([zm-2,zm+2]);
%     zticks = [round(zm)-2:1:round(zm)+2];
%     set(gca,'ZTick',zticks);
%     set(gca,'ZTickLabel',sprintf('%+3.0f\n',zticks))
%     tightInset = get(gca, 'TightInset');
%     position(1) = tightInset(1);
%     position(2) = tightInset(2);
%     position(3) = 1 - tightInset(1) - tightInset(3);
%     position(4) = 1 - tightInset(2) - tightInset(4);
%     set(gca, 'Position', position);
%     set(f,'Units','Inches');
%     pos = get(f,'Position');
%     set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    
    %Options for Sphere Simulation
    axis image;
    zcs = sort(Zp(2,:));
    zmean = mean(zcs(1:ceil(length(zcs)/2)));
    axis([-50 50 -50 50 zmean-(50*2) zmean+(50*2)]);
    
%     zmean = min(Zp(2,:));
%     axis([-8 8 -8 8 zmean-4 zmean+16]);
    
    zticks = get(gca,'ZTick');
    set(gca,'ZTickLabel', sprintf('%+5.0f\n',zticks))
    
    tightInset = get(gca, 'TightInset');
    position(1) = tightInset(1);
    position(2) = tightInset(2);
    position(3) = 1 - tightInset(1) - tightInset(3);
    position(4) = 1 - tightInset(2) - tightInset(4);
    set(gca, 'Position', position);
    set(f,'Units','Inches');
    pos = get(f,'Position');
    set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    
    saveas(f, output);
    
    index = index + 1
    
    %close all;
end
