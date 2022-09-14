% smooth functions
% 输入 z_ik  p_ik(用户定义误差)  ni(智能体阶次)  k(智能体第几阶次的状态量)
% 输出 sg_ik(z_ik)

function sg_function = sg_ik(z_ik,p_ik,ni,k)
    if abs(z_ik) >= p_ik
        sg_function = z_ik/abs(z_ik);
    else
        sg_function = z_ik/(((p_ik^2 - z_ik^2)^(ni-k+2)) + abs(z_ik));
    end
end

