% f_ik
% ÊäÈë z_ik p_ik
% Êä³ö f_ik(z_ik)

function f_ik_function = f_ik(z_ik,p_ik)
    if abs(z_ik) >= p_ik
        f_ik_function = 1;
    else 
        f_ik_function = 0;
    end
end

