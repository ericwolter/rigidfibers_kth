input=['XcT_gen1000.in']; 
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

N = size(pos,1);

domain = [min(min(Xp)),min(min(Yp)),min(min(Zp));max(max(Xp)),max(max(Yp)),max(max(Zp))]

nearest_distance_segment = 99999999.0;
nearest_distance_center = 99999999.0;
minimal_distance_segment = 99999999.0;
minimal_distance_center = 99999999.0;
distance_segment = 0.0;
distance_center = 0.0;
total_segment = 0.0;
total_center = 0.0;

for i = 1:N
    P1 = [Xp(1,i) Yp(1,i) Zp(1,i)];
    C = [Xp(2,i) Yp(2,i) Zp(2,i)];
    P2 = [Xp(3,i) Yp(3,i) Zp(3,i)];
    for j = 1:N
        if i ~= j
            P3 = [Xp(1,j) Yp(1,j) Zp(1,j)];
            D = [Xp(2,j) Yp(2,j) Zp(2,j)];
            P4 = [Xp(3,j) Yp(3,j) Zp(3,j)];
            
            distance_segment = DistBetween2Segment(P1,P2,P3,P4);
            distance_center = norm(D-C);
            
            if distance_segment < nearest_distance_segment
                nearest_distance_segment = distance_segment;
            end
            if distance_center < nearest_distance_center
                nearest_distance_center = distance_center;
            end
        end
    end
    
    total_segment = total_segment + nearest_distance_segment;
    total_center = total_center + nearest_distance_center;
    
    if nearest_distance_segment < minimal_distance_segment
        minimal_distance_segment = nearest_distance_segment;
    end
    if nearest_distance_center < minimal_distance_center
        minimal_distance_center = nearest_distance_center;
    end
end

minimal_distance_segment = minimal_distance_segment
minimal_distance_center = minimal_distance_center
average_distance_segment = total_segment / N
average_distance_center = total_center / N

    
