%%% SLAM simulator

%% initilaize the simulator
init_simulator();
landmarks = importdata('map.txt',' ');
% initialize values
mu = [0,0,0]';
mu_all = mu(1:2,:)';
sigma = diag([0.01, 0.01, 0.001]);
N = [];
u = [0;0];
unchanged = 0;
limit_xy = [-37, 37];
turn_cnt = 1;
max_cnt = 1;
bool_cnt = true;

%% start the simulation
while true
    %% get the direction of the robot  
  turn_cnt = turn_cnt - 1;
  u = get_direction(mu, turn_cnt);
  if turn_cnt == 0
    if bool_cnt == true
      bool_cnt = false;
      turn_cnt = max_cnt;
    else
      max_cnt = max_cnt + 1;
      turn_cnt = max_cnt;
      bool_cnt = true;
    end
  end
  
  %% translate the robot
  move_robot(u(1), 0);
  % and get landmark
  max = 0;
  l = [];
  for i=1:10
    l_tmp = get_landmarks();
    if size(l_tmp,2) > max
      l = l_tmp;
      max = size(l_tmp,2);
    end
  end
  
  %% perform EKF_SLAM_UC
  u_t = [u(1); 0];
  [mu, sigma, N] = ekf_slam_uc(mu, sigma, u_t, l);
  mu_all = [mu_all; mu(1:2,1)'];
  
  %% rotate the robot
  move_robot(0, u(2));
  % and get landmark
  max = 0;
  l = [];
  for i=1:10
    l_tmp = get_landmarks();
    if size(l_tmp,2) > max
      l = l_tmp;
      max = size(l_tmp,2);
    end
  end

  %% perform EKF_SLAM_UC
%   N_prev = N;
  u_r = [0; u(2)];
  [mu, sigma, N] = ekf_slam_uc(mu, sigma, u_r, l);
  mu_all = [mu_all; mu(1:2,1)'];
  
  %% visualize
  figure(1); clf;
  hold on;
  t=-pi:0.01:pi;
  x_mu = mu(1) + sqrt(2*sigma(1,1))*cos(t);
  y_mu = mu(2) + sqrt(2*sigma(2,2))*sin(t);
  plot(x_mu, y_mu, 'r');
  for i = 4:2:size(mu,1)
    x = mu(i) + sqrt(2*sigma(i,i))*cos(t);
    y = mu(i+1) + sqrt(2*sigma(i+1,i+1))*sin(t);
    plot(x, y, 'g');
  end
  plot(mu(1), mu(2), 'rs');
  plot(mu_all(:,1), mu_all(:,2), 'r-', mu(4:2:end), mu(5:2:end), 'g*', landmarks(:,1), landmarks(:,2), 'bo');
  xlim([-60 60]); ylim([-60 60]);
  hold off;
  pause(1);
  
  % check boundary
%   if mu(1) < limit_xy(1)*2 || mu(2) < limit_xy(1)*2 || mu(1) > limit_xy(2)*2 || mu(2) > limit_xy(2)*2
%     break;
%   end
  % check error
%   if sqrt(sigma(1,1)) > 5 || sqrt(sigma(2,2)) > 5 || sqrt(sigma(3,3)) > pi/16
%     break;
%   end
%   % if N is not changed for the certain amount of times, just break;
%   if N == N_prev
%     unchanged = unchanged + 1;
%   else
%     unchanged = 0;
%   end
%   
%   if unchanged > 10
%     break;
%   end
  
end % end of the simulation
