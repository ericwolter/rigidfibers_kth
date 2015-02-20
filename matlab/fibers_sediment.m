setup = 'sphere_2000_0040';
files = dir(['input_' setup '/*.out']);
% files = files(1455:end);
close all;
index = 0;
% f = figure();
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperUnits', 'inches');
% set(gcf, 'PaperSize', [17 22]);
% set(gcf, 'PaperPosition', [1 1 15 20]);

input=['input_' setup '/' files(1).name];
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

N_0 = size(Xp,2);

% Calculate center of mass
FC = [Xp(2,:);Yp(2,:);Zp(2,:)];
CENTER_OF_MASS = mean(FC,2);

FC_prev = zeros(size(FC));
CENTER_OF_MASS_PREV = zeros(size(CENTER_OF_MASS));

HORIZONTAL_DIST = sqrt(sum((FC(1:2,:) - repmat(CENTER_OF_MASS(1:2),1,N_0)).^2));
R_0 = max(HORIZONTAL_DIST);
blacklist = (HORIZONTAL_DIST > R_0);
% Xp(:,blacklist)=[];
% Yp(:,blacklist)=[];
% Zp(:,blacklist)=[];

% size(Xp)

% plot3(Xp,Yp,Zp,'LineWidth',0.5);
% grid on;
% axis image;

PLOT_VX = [];
PLOT_VY = [];
PLOT_VZ = [];

PLOT_RX = [];
PLOT_RY = [];
PLOT_RZ = [];
PLOT_SX = [];
PLOT_SY = [];
PLOT_SZ = [];
PLOT_N = [];

for file = files'
    input=['input_' setup '/' file.name]; 
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
    
    Xg = Xp;
    Yg = Yp;
    Zg = Zp;
    
    Xg(:,blacklist)=[];
    Yg(:,blacklist)=[];
    Zg(:,blacklist)=[];
    
    % Calculate center of mass
    N = size(Xg,2);
    PLOT_N = [PLOT_N N/N_0];
    
    ALL_CENTERS = [Xp(2,:);Yp(2,:);Zp(2,:)];
    
    ALL_VELOCITIES = ALL_CENTERS - FC_prev;
    FC_prev = ALL_CENTERS;
    ALL_VELOCITIES(:,blacklist)=[];
    
    PLOT_VY = [PLOT_VY std(ALL_VELOCITIES(3,:))];
    PLOT_VZ = [PLOT_VZ mean(ALL_VELOCITIES(3,:))/(N_0/R_0)];
    
    CLOUD_CENTERS = [Xg(2,:);Yg(2,:);Zg(2,:)];
    CENTER_OF_MASS = mean(CLOUD_CENTERS,2);
    
    CENTER_VELOCITY = CENTER_OF_MASS - CENTER_OF_MASS_PREV;
    CENTER_OF_MASS_PREV = CENTER_OF_MASS;
    PLOT_VX = [PLOT_VX CENTER_VELOCITY(3)/(N_0/R_0)];
    
    DIFF_CLOUD = CLOUD_CENTERS - repmat(CENTER_OF_MASS,1,N);
    HORIZONTAL_DIST_CLOUD = sqrt(sum((DIFF_CLOUD(1:2)).^2));
    R = max(HORIZONTAL_DIST_CLOUD);
    
    PLOT_RX = [PLOT_RX max(abs(DIFF_CLOUD(1,:)))/R_0];
    PLOT_RY = [PLOT_RY max(abs(DIFF_CLOUD(2,:)))/R_0];
    PLOT_RZ = [PLOT_RZ max(abs(DIFF_CLOUD(3,:)))/R_0];
    PLOT_SX = [PLOT_SX std(abs(DIFF_CLOUD(1,:))/R_0)];
    PLOT_SY = [PLOT_SY std(abs(DIFF_CLOUD(2,:))/R_0)];
    PLOT_SZ = [PLOT_SZ std(abs(DIFF_CLOUD(3,:))/R_0)];
    
    VERTICAL_DIST_ALL = abs(ALL_CENTERS(3,:) - repmat(CENTER_OF_MASS(3),1,N_0));
    blacklist = (VERTICAL_DIST_ALL > R_0);
    
%     zcs = sort(Zp(2,:));
%     zinterest = zcs(1:ceil(length(zcs)/2));
%     zmean = mean(zinterest);
%     zstd = std(zinterest);

%     subplot(2,2,1);
%     plot3(Xp,Yp,Zp,'LineWidth',0.5);
%     title('Orthographic');
%     grid on;
%     axis image;
%     axis([-8 8 -8 8 zmean-(8*2) zmean+(8*2)]);
%     
%     subplot(2,2,2);
%     plot3(Xp,Yp,Zp,'LineWidth',0.5);
%     title('Front');
%     %x
%     az = 90;
%     el = 0;
%     view(az, el);
%     grid on;
%     axis image;
%     axis([-8 8 -8 8 zmean-(8*2) zmean+(8*2)]);
%     
%     subplot(2,2,3);
%     plot3(Xp,Yp,Zp,'LineWidth',0.5);
%     title('Left');
%     % y
%     az = -90;
%     el = 0;
%     view(az, el);
%     grid on;
%     axis image;
%     axis([-8 8 -8 8 zmean-(8*2) zmean+(8*2)]);
% 
%     subplot(2,2,4);
%     plot3(Xp,Yp,Zp,'LineWidth',0.5);
%     title('Top');
%     % z
%     az = 1e-10;
%     el = 90;
%     view(az, el);
%     grid on;
%     axis image;
%     axis([-8 8 -8 8 zmean-(8*2) zmean+(8*2)]);
    
%     f = figure(1);
%     plot3(Xp,Yp,Zp,'LineWidth',0.5);
%     grid on;
%     axis image;
%     axis([-10 10 -10 10 CENTER_OF_MASS(3)-(10*2) CENTER_OF_MASS(3)+(10*2)]);
%     output=['output_' setup '/state_' sprintf('%05d',index) '.pdf'];
    
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
%     axis image;
%     zcs = sort(Zp(2,:));
%     zmean = mean(zcs(1:ceil(length(zcs)/2)));
%     axis([-50 50 -50 50 zmean-(50*2) zmean+(50*2)]);
    
%     zmean = min(Zp(2,:));
%     axis([-8 8 -8 8 zmean-4 zmean+16]);
    
%     zticks = get(gca,'ZTick');
%     set(gca,'ZTickLabel', sprintf('%+5.0f\n',zticks))
    
%     tightInset = get(gca, 'TightInset');
%     position(1) = tightInset(1);
%     position(2) = tightInset(2);
%     position(3) = 1 - tightInset(1) - tightInset(3);
%     position(4) = 1 - tightInset(2) - tightInset(4);
%     set(gca, 'Position', position);
%     set(f,'Units','Inches');
%     pos = get(f,'Position');
%     set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    
%      saveas(f, output);
    
    index = index + 1
    
    %close all;
end

figure(1);
plot(smooth(PLOT_RX,10));
hold on;
plot(smooth(PLOT_RY,10));
hold on;
plot(smooth(PLOT_RZ,10));
title('Radius Max');

figure(2);
plot(PLOT_SX);
hold on;
plot(PLOT_SY);
hold on;
plot(PLOT_SZ);
title('Radius Std');

figure(3);
plot(PLOT_SZ);
hold on;
plot(smooth(PLOT_SZ,20));
hold on;
plot(smooth(PLOT_SZ,100));
hold on;
title('R_{SZ} Std');

figure(4);
plot(PLOT_N);
hold on;
title('Number of fibers');

figure(5);
plot(PLOT_VX);
hold on;
plot(PLOT_VZ);
title('Velocity');

% csvwrite('sphere_2000_0040.csv',[[1:1:size(PLOT_N,2)]' PLOT_N' PLOT_SZ' PLOT_VZ'])

% f = figure(3);
% plot(diff(smooth(diff(smooth(PLOT_RX,10)),5)));
% hold on;
% plot(diff(smooth(diff(smooth(PLOT_RY,10)),5)));
% title('Radius Accel');
% 
% f = figure(4);
% plot(PLOT_N);
% title('Number');
% f = figure(5);
% plot(diff(diff(PLOT_N)));
% title('Number Accel');

if find(PLOT_N<0.97,1) > 1
    start_index = find(PLOT_N<0.97,1);
else
    start_index = 1;
end

if find(PLOT_N<0.5,1) > 1
    end_index = find(PLOT_N<0.5,1);
else
    end_index = index;
end
[~, min_index] = min(PLOT_SZ(start_index:end_index));
[~, max_index] = max(PLOT_SZ(start_index:end_index));
min_index = start_index + min_index - 1
max_index = start_index + max_index - 1
mid_index = round((max_index + min_index) / 2)

% f = figure(3);
% plot(PLOT_N);

