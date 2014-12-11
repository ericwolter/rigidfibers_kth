files = dir('input/*.out');

index = 0;

f = figure();
for file = files'
    input=['input/' file.name]; 
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
    
    plot3(Xp,Yp,Zp);
    grid on;
    axis([-30 30 -30 30 -1000 10]);
    
    output=['output/' sprintf('%05d',index) '.png'];
    
    set(f,'PaperUnits','inches','PaperPosition',[0 0 5 8])
    saveas(f, output);
    
    index = index + 1
    
    %close all;
end

