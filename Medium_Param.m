%% Medium parameter
density = 998;
p0 = 101.325 * 10^3; % equallibrium pressure 101.325 kPa (standard atmosphere pressure)
water_surf_tens = 72*10^-3;  %0.072 N/m
c = 1480; % speed of sound
k = 1.095; % gas polytropic index
% p_v = 2.266 * 10^3; % vapour pressure 226 kPa
p_v = 0;
viscosity = 0.001; % temperature at 20 celcius (shear viscosity of surrounding liquid)