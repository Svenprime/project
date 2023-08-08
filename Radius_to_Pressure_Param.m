%% R''
run('Medium_Param.m');
run('Bubble_Shell_Param.m');
R_buckling = 0.91 * R0;
R_break = 1.08 * R0;
equ_surf_tens = elas_mod *(R0^2 / R_buckling^2 - 1);%effective surface tension constant value
for i = 1:length(Rt)

    if Rt(i) < R_buckling
        surf_tens = 0;
    elseif Rt(i) <= R_break
        surf_tens = elas_mod *(Rt(i)^2 / R_buckling^2 - 1);
    else
        surf_tens = water_surf_tens;
    end

    %the acceleration derived from equation3 in the RPM equation
    R_a(i,1) = (-3*bub_veloc(i)^2 / 2 + 1/density ...
        * ((p0-p_v + 2*equ_surf_tens / R0) * (R0/Rt(i))^(3*k) ...%equation 7 the internal pressure
        *(1-(3*k* bub_veloc(i)/c)) - ((2*surf_tens) / Rt(i)) ...
        + p_v - 4*viscosity* bub_veloc(i)/Rt(i) - 4*k_s* bub_veloc(i)/Rt(i)^2 - p0 - p_r(i)')) / Rt(i);
    
end