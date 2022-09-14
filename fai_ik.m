% RBFNN�����磬���Աƽ������Ժ���
% ��_ik(��_ik)
% ����Ϊ niΪ������״̬�� kΪ�ڼ��״�״̬ �Ա�������x(����_i,j)
% ���Ϊ 

% norm Ϊ����
% aΪ�˺������� 
function fai_ik = fai_ik(x,ni,k,center)                
    widths =  8^2;                             % ��˹�˺����Ļ���
    a = center * ones(size(x)) ;
    fai_ik = zeros(1,2*ni);
    fai_ik(k) = exp( - (norm(x-a))^2 / widths );
    fai_ik(ni+k) = 1;
end

