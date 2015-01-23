close all;

setup = 'force2_2000_0040_2';
files = dir(['input_' setup '/*.out']);
index = 0;

co_reds = [
    0.99608,0.89804,0.85098;
    0.98824,0.73333,0.63137;
    0.98824,0.57255,0.44706;
    0.98431,0.41569,0.2902;
    0.93725,0.23137,0.17255;
    0.79608,0.09412,0.11373;
    0.6,0,0.05098
];

co_blues = [
    0.93725,0.95294,1;
    0.77647,0.85882,0.93725;
    0.61961,0.79216,0.88235;
    0.41961,0.68235,0.83922;
    0.25882,0.57255,0.77647;
    0.12941,0.44314,0.7098;
    0.03137,0.27059,0.58039
];


f = figure();

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
co_reds = repmat(co_reds,N_0,1);
co_blues = repmat(co_blues,N_0,1);

% Calculate center of mass
FC = [Xp(2,:);Yp(2,:);Zp(2,:)];
CENTER_OF_MASS = mean(FC,2);

FC_prev = zeros(size(FC));
CENTER_OF_MASS_PREV = zeros(size(CENTER_OF_MASS));

HORIZONTAL_DIST = sqrt(sum((FC(1:2,:) - repmat(CENTER_OF_MASS(1:2),1,N_0)).^2));
R_0 = max(HORIZONTAL_DIST);
blacklist = (HORIZONTAL_DIST > R_0);

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
    
    h = plot3(Xp,Yp,Zp,'LineWidth',0.5);
    for m=1:N_0/2
        set(h(m),'Color',co_reds(m,:))
    end
    for m=N_0/2:N_0
        set(h(m),'Color',co_blues(m,:))
    end
    
    grid on;
    axis 'image';
    axis([-R_0*3 R_0*3 -R_0*3 R_0*3 CENTER_OF_MASS(3)-(R_0*2) CENTER_OF_MASS(3)+(R_0*2*3)]);
    
%     subplot(2,2,2);
%     plot3(Xp,Yp,Zp,'LineWidth',0.5);
%     title('Front');
%     %x
%     az = 90;
%     el = 0;
%     view(az, el);
%     grid on;
%     axis image;
%     axis([-R_0*2 R_0*2 -R_0*2 R_0*2 CENTER_OF_MASS(3)-(R_0*2*2) CENTER_OF_MASS(3)+(R_0*2*2)]);
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
%     axis([-R_0*2 R_0*2 -R_0*2 R_0*2 CENTER_OF_MASS(3)-(R_0*2*2) CENTER_OF_MASS(3)+(R_0*2*2)]);
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
%     axis([-R_0*2 R_0*2 -R_0*2 R_0*2 CENTER_OF_MASS(3)-(R_0*2*2) CENTER_OF_MASS(3)+(R_0*2*2)]);
    
    zticks = get(gca,'ZTick');
    set(gca,'ZTickLabel', sprintf('%+5.0f\n',zticks))
    
%     tightInset = get(gca, 'TightInset');
%     position(1) = tightInset(1);
%     position(2) = tightInset(2);
%     position(3) = 1 - tightInset(1) - tightInset(3);
%     position(4) = 1 - tightInset(2) - tightInset(4);
%     set(gca, 'Position', position);
%     pos = get(f,'Position');
%     set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    
    output=['output_' setup '/state_' sprintf('%05d',index) '.pdf'];
    saveas(f, output);
    
    index = index + 1
    
    %close all;
end

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
