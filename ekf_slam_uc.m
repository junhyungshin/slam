function [mu, sigma, N] = ekf_slam_uc(mu, sigma, u, l)
  % constants
  deltaT = 1;
  alpha = 28;
  big_val = 5764801; % 7^8

  z = [];
  for i=1:size(l,2)
    z(1,i) = l(i).d;
    z(2,i) = wrapToPi(l(i).b);
  end

  % prediction
  N = (size(mu,1)-3)/2;
  F = [eye(3), zeros(3,2*N)];
  if u(2) == 0
%     mu_bar = mu + F'*[u(1)*deltaT*cos(mu(3));
%                       u(1)*deltaT*sin(mu(3));
%                       u(2)*deltaT];
    mu_bar = mu + F'*[u(1)*deltaT*cos(mu(3));
                      u(1)*deltaT*sin(mu(3));
                      0];
    G = eye(3+2*N) + F'*[0,0,-1*u(1)*deltaT*sin(mu(3));
                         0,0,u(1)*deltaT*cos(mu(3));
                         0,0,0]*F;                    
  elseif u(1) == 0
%       mu_bar = mu + F'*[-1*u(1)/u(2)*sin(mu(3))+u(1)/u(2)*sin(mu(3)+u(2)*deltaT);
%                         u(1)/u(2)*cos(mu(3))-u(1)/u(2)*cos(mu(3)+u(2)*deltaT);
%                         u(2)*deltaT];
    mu_bar = mu + F'*[0; 0; u(2)*deltaT];
    mu_bar(3) = wrapToPi(mu_bar(3));
%     G = eye(3+2*N) + F'*[0,0,u(1)/u(2)*cos(mu(3))-u(1)/u(2)*cos(mu(3)+u(2)*deltaT);
%                          0,0,u(1)/u(2)*sin(mu(3))-u(1)/u(2)*sin(mu(3)+u(2)*deltaT);
%                          0,0,0]*F;
    G = eye(3+2*N) + F'*zeros(3,3)*F;
  end

  R = diag([0.001, 0.001, 0.001]);
  sigma_bar = G*sigma*G' + F'*R*F;
  Q = diag([0.01, 0.001]);

  % correction
  for i = 1:size(z,2)
    mu_bar_obs = [mu_bar(1); mu_bar(2)] + z(1,i)*[cos(z(2,i) + mu_bar(3)); sin(z(2,i) + mu_bar(3))];
    mu_bar = [mu_bar; mu_bar_obs];
    sigma_bar = [[sigma_bar, zeros(3+2*N,2)]; [zeros(2,3+2*N), diag([big_val, big_val])]];
    for k = 4:2:size(mu_bar, 1)
      k_tmp = (k-4)/2+1; % 1,2,3...
      delta_k = [mu_bar(k)-mu_bar(1); mu_bar(k+1)-mu_bar(2)];
      q_k = delta_k'*delta_k;
      z_hat(:,k_tmp) = [sqrt(q_k); wrapToPi(atan2(delta_k(2), delta_k(1))-mu_bar(3))];
      F_k = [[eye(3), zeros(3,2*(k_tmp)-2), zeros(3,2), zeros(3,2*(N+1)-2*(k_tmp))];
             [zeros(2,3), zeros(2,2*(k_tmp)-2), eye(2), zeros(2,2*(N+1)-2*(k_tmp))]];
      H{k_tmp} = 1/q_k*[-1*sqrt(q_k)*delta_k(1), -1*sqrt(q_k)*delta_k(2), 0, sqrt(q_k)*delta_k(1), sqrt(q_k)*delta_k(2);
                    delta_k(2), -1*delta_k(1), -q_k, -1*delta_k(2), delta_k(1)]*F_k;
      psi(:,:,k_tmp) = H{k_tmp}*sigma_bar*(H{k_tmp}')+Q;
      z_tmp = z(:,i)-z_hat(:,k_tmp);
      z_tmp(2) = wrapToPi(z_tmp(2));
      pi(k_tmp) = (z_tmp)'*inv(psi(:,:,k_tmp))*(z_tmp);      
    end
    pi(N+1) = alpha;
    j(i) = find(pi == min(pi));
    N_prev = N;
    N = max(N, j(i));
    if N == N_prev
      mu_bar = mu_bar(1:3+2*N);
      sigma_bar = sigma_bar(1:3+2*N,1:3+2*N);      
      H{j(i)} = H{j(i)}(:,1:end-2);
    end
    K{i} = sigma_bar*H{j(i)}'*inv(psi(:,:,j(i)));
    z_tmp = z(:,i)-z_hat(:,j(i));
    z_tmp(2) = wrapToPi(z_tmp(2));
    mu_bar = mu_bar + K{i}*(z_tmp);
    kh = K{i}*H{j(i)};
    sigma_bar = (eye(size(kh)) - kh)*sigma_bar;
  end
  mu = mu_bar;
  sigma = sigma_bar;
  
  
  
end