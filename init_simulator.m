function init_simulator()

% (re)initializes the simulator

clear move_robot % removes persistent variables
clear get_landmarks

x= 0;
y =0;
theta = 0;

save('swap.mat','x','y','theta')