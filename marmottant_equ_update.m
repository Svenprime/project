function drdt = marmottant_equ_update(t,R,f,cycle,d_p,T,phase,R0)
%% Applied Acoustic Force
if t < T*cycle
    Rasp = d_p*sin(2*pi*f*t+phase);
else
    Rasp =0;
end
p_r = Rasp;
%% Liquid Parameter
run('Medium_Param.m');
%%
run('Bubble_Shell_Param.m');
%% Marmottant Model surface tension
R_buckling = 0.91 * R0;
R_break = 1.08 * R0;
equ_surf_tens = elas_mod *(R0^2 / R_buckling^2 - 1);
if R(1) <= R_buckling
    surf_tens = 0;
elseif R(1) <= R_break
    surf_tens = elas_mod *(R(1)^2 / R_buckling^2 - 1);
else
    surf_tens = water_surf_tens;
end

%% Solve Equation        
RP_Equ = (p0 + 2 * equ_surf_tens / R0 - p_v) * (R0/R(1))^(3*k)*(1-(3*k*R(2)/c)) - (p0-p_v + p_r);
Surf_Tens_Damping = 2*surf_tens / R(1);
Visco_Damping = 4*viscosity*R(2)/R(1) + 4*k_s*R(2)/R(1)^2;
RP_Equ_Marmot = ((RP_Equ - Surf_Tens_Damping - Visco_Damping) / density - (3*R(2)^2 / 2)) / R(1);
drdt = [R(2); RP_Equ_Marmot];
end