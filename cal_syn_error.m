% ����������ͬ�����  z_i,1
% ����AΪ�ڽӾ���
%     BΪ��leader�����Ӿ���
%     xΪ�������ʼ״̬
%     ydΪleader״̬
%     follower_numΪ����������
% ���Ϊ ei����
function E = cal_syn_error(A,B,x,yd,follower_num)  
    E = zeros(follower_num,1);
    E_ij = zeros(follower_num,follower_num);
    for i = 1 : follower_num
        for j = 1 : follower_num
            E_ij(i,j) = A(i,j)*(x(i)-x(j));
        end
    end
    E_adjacency = sum(E_ij,2);     % ���ڽ�����������(�������)

    for i = 1 : follower_num
        E(i) = B(i)*(x(i)-yd) + E_adjacency(i);  % leader�����+�ڽ����������
    end
end

