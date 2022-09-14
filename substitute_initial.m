% ����ʼ״̬��������һʱ�̵�״̬����
%   ���룺x Ϊ��������ʱ�̵�״̬����
%         y Ϊleader��״̬
%         a Ϊÿ��������״̬�״Σ���ÿ���������ж��ٸ�״̬����
%         b Ϊ������ĸ���
%         t ���沽��

function first_time_state = substitute_initial(x,y,a,b,t)
    a = a + 4;                          % ״̬������Ҫ����leader������͵���
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

