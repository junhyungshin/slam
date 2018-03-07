
function move_robot(ut,ur)

% moves the robot applying translational speed ut and rotational speed ur
% for 1 second. Only one of the two can be different from 0, i.e., either
% you move forward or you rotate, but you do not move along an arc

persistent x y theta

% initialize persistent variables when needed
if isempty(x)
    x =0;
    y= 0;
    theta = 0;
end

deltaT= 1;

varx = 0.01;
vart = 0.001;

x = x + ut * deltaT * cos(theta) + sqrt(varx)*randn(1);
y = y + ut * deltaT * sin(theta)+ sqrt(varx)*randn(1);
theta = theta + ur*deltaT+ sqrt(vart)*randn(1);

save('swap.mat','x','y','theta')

end
