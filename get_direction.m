%% get direction of the robot by return u_t and u_r
function u = get_direction(mu, turn_cnt)
  move_unit = 7; % robot always moves 7m
  
  if turn_cnt == 0
    u_r = pi/2;
  else
    u_r = 0;
  end
  u_t = move_unit;
  u = [u_t, u_r];
end