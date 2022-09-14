% RBFNN神经网络，用以逼近非线性函数
% φ_ik(χ_ik)
% 输入为 ni为智能体状态数 k为第几阶次状态 自变量向量x(即χ_i,j)
% 输出为 

% norm 为求范数
% a为核函数中心 
function fai_ik = fai_ik(x,ni,k,center)                
    widths =  8^2;                             % 高斯核函数的基宽
    a = center * ones(size(x)) ;
    fai_ik = zeros(1,2*ni);
    fai_ik(k) = exp( - (norm(x-a))^2 / widths );
    fai_ik(ni+k) = 1;
end

