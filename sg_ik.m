% smooth functions
% ���� z_ik  p_ik(�û��������)  ni(������״�)  k(������ڼ��״ε�״̬��)
% ��� sg_ik(z_ik)

function sg_function = sg_ik(z_ik,p_ik,ni,k)
    if abs(z_ik) >= p_ik
        sg_function = z_ik/abs(z_ik);
    else
        sg_function = z_ik/(((p_ik^2 - z_ik^2)^(ni-k+2)) + abs(z_ik));
    end
end

