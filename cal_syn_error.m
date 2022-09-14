% 计算智能体同步误差  z_i,1
% 输入A为邻接矩阵
%     B为与leader的连接矩阵
%     x为智能体初始状态
%     yd为leader状态
%     follower_num为智能体数量
% 输出为 ei矩阵
function E = cal_syn_error(A,B,x,yd,follower_num)  
    E = zeros(follower_num,1);
    E_ij = zeros(follower_num,follower_num);
    for i = 1 : follower_num
        for j = 1 : follower_num
            E_ij(i,j) = A(i,j)*(x(i)-x(j));
        end
    end
    E_adjacency = sum(E_ij,2);     % 与邻接智能体的误差(按行求和)

    for i = 1 : follower_num
        E(i) = B(i)*(x(i)-yd) + E_adjacency(i);  % leader的误差+邻接智能体误差
    end
end

