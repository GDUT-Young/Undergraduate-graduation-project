function D = cal_indegree_matrix(A)
% �����ڽӾ��������Ⱦ���
    [row,column] = size(A);  %�����ڽӾ�������row ����column
    d = zeros(1,row);
    sum_A = sum(A,2);      % �ڽӾ���ÿһ�����
    for i = 1 : row
        d(i) = sum_A(i);
    end
    D = diag(d,0);       % DΪ��Ⱦ���
end

