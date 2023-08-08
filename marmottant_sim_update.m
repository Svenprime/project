clearvars
%% Acoustic Parameter
run('Acoustic_Param.m');


%% Rayleigh-Plesset-Marmottant model ode45 solver

tspan = [0 5*10^-6];
[t,R] = ode45(@(t,R) marmottant_equ_update(t,R,f,cycle,d_p,T,phase,R0),tspan,[R0;0]); % time 1 mm

Rt = R(:,1); % Radius
bub_veloc = R(:,2);  % Bubble Wall Veloc


%% Plot Bubble Radius Curve

Rasp = d_p*sin(2*pi*f*t'+phase); % Applied Acoustic Force
for i = 1 : length(t)
    
    if t(i,1) > T * cycle
       Rasp(1,i) = 0;
    end
end


figure;

subplot(311)
plot(t*10^6,(Rasp)*10^-3)
xlabel('time(\mus)');
ylabel('Pressure(kPa)');
title('Driving Pulse')

subplot(312)
plot(t*10^6,Rt*10^6);
xlabel('time(\mus)');
ylabel('Radius(\mum)');
title('Bubble Radius')

subplot(313)
plot(t*10^6,bub_veloc);
title('Bubble Wall Velocity')
xlabel('time(\mus)');
ylabel('velocity(m/s)');

sgtitle(sprintf('Rayleigh-Plesset Equation Bubble Response in %d kPa', d_p/10^3))

%% Radius to Pressure
p_r = Rasp;
run('Radius_to_Pressure_Param.m');
r = 3*10^-3; % 3mm distance from the center of cavity
sw_p = density * (Rt.^2.*R_a + 2*Rt.*bub_veloc.^2)/r; % Acoustic Emission Calculation equation8
%% Plot Derived Pressure Waveform  
figure;
plot(t,sw_p*10^-3)
xlabel('time(\mus)');
ylabel('Pressure(kPa)');
title('Shockwave waveform');
title(sprintf('Shockwave waveform in %d kPa', d_p/10^3))
