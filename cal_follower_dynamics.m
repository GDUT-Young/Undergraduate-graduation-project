% follower's dynamics
% 计算智能体的状态
% 输入： x 迭代前状态
%        u 输入
%        t 仿真时间
% 输出： follower_dynamics 加入控制率后的智能体状态
% xi_k_j+1 = xi_k_j +  ( xi_k_j_dot * dt )    状态下一时刻的值等于上一时刻值加上其

function follower_dynamics = cal_follower_dynamics(x,u,t,dt)
    follower_dynamics = zeros(1,6);  % the column represent xi_1 xi_2 xi_1_dot xi_2_dot yd yd_dot

    follower_dynamics(1) = x(1) + ( x(3) * dt );                     % xi_1 = yi 智能体输出
    follower_dynamics(2) = x(2) + ( x(4) * dt );                     % xi_2      
    
    follower_dynamics(4) = u + follower_dynamics(2) *...
                                   sin(0.2/(4 + follower_dynamics(1)^2));             % xi_2_dot    
    follower_dynamics(3) = follower_dynamics(2) +...
                        ( (follower_dynamics(2))^3 / (2 + ( follower_dynamics(1))^2 ) );  % xi_1_dot
 
    follower_dynamics(5) = x(5) + ( x(6) * dt );                                      % yd
    follower_dynamics(6) = -1.4*follower_dynamics(5) + 2.5*cos(0.6*t)+2*(t^2)*exp(-t)-0.5; % yd_dot  
end

