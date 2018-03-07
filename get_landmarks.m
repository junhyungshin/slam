function l = get_landmarks()

persistent list

if isempty(list)
    list = importdata('map.txt',' ');
end

cutoff = 7;

coords = load('swap.mat');
x = coords.x;
y = coords.y;
theta = coords.theta;

rangevar = 0.01;
bearvar = 0.001;
missed_detection_prob = 0.1;

l = [];

for i =1:size(list,1)
   
    d = hypot(list(i,1)-x,list(i,2)-y);
    if d < cutoff
        if rand(1) > missed_detection_prob % random missed detection
            reading.d = d + sqrt(rangevar)*randn(1);
            u = [cos(theta) sin(theta) 0];
            v = [list(i,1)-x list(i,2)-y 0];
            g = cross(u,v);
            if g(3)>0
                reading.b = atan2(norm(cross(u,v)),dot(u,v));
            else
                reading.b = -atan2(norm(cross(u,v)),dot(u,v));
            end
            reading.b = reading.b+ sqrt(bearvar)*randn(1);
            l = [ l reading];
        end
    end
end

