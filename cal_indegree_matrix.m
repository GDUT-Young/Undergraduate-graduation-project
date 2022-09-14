function D = cal_indegree_matrix(A)
% 根据邻接矩阵计算入度矩阵
    [row,column] = size(A);  %计算邻接矩阵行数row 列数column
    d = zeros(1,row);
    sum_A = sum(A,2);      % 邻接矩阵每一行求和
    for i = 1 : row
        d(i) = sum_A(i);
    end
    D = diag(d,0);       % D为入度矩阵
end

