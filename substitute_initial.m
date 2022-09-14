% 将初始状态代入计算第一时刻的状态矩阵
%   输入：x 为智能体零时刻的状态矩阵
%         y 为leader的状态
%         a 为每个智能体状态阶次（即每个智能体有多少个状态量）
%         b 为智能体的个数
%         t 仿真步长

function first_time_state = substitute_initial(x,y,a,b,t)
    a = a + 4;                          % 状态的数量要加上leader的输出和导数
    first_time_state = zeros(b,a);
    for i = 1 : b        
        first_time_state(i,1) = x(i,1);      % xi_1
        first_time_state(i,2) = x(i,2);      % xi_2
        first_time_state(i,5) = y;           % yd
        
        first_time_state(i,3) = first_time_state(i,2) + ((first_time_state(i,2)^3)/...
                                    (2 + (first_time_state(i,1)^2)));     % xi_1_dot
        first_time_state(i,4) = 0 + first_time_state(i,2) * sin(0.2/(4 + (first_time_state(i,1)^2))); % xi_2_dot
        first_time_state(i,6) = -1.4*first_time_state(i,5) + 2.5*cos(0.6*t)+2*(t^2)*exp(-t)-0.5;   % yd_dot
    end
end

