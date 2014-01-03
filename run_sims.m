



time = 1:1/60:40;                                    % Time steps for simulation [hrs]
setpoint = [time', 70*ones(length(time),1)];        % Set point for normal operations [deg F]
outdoor_temp = [time', 5*sin(2*pi/24*time') + 60]; % Outdoor temperature [deg F]

learning_mode = 'On';   % Status of the learning algorithm ('On' or 'Off')


heater_sf = 0.8:0.05:1.2;

% Heater temperature sensitivity
c_error_heat = zeros(length(heater_sf),4);

for i = 1:length(heater_sf)
    [t, heater_status, room_temp, c_actual, c_estimated] = thermSim(time, setpoint, outdoor_temp, 'On',0, 0, 1, heater_sf(i));
    t_heater_estimated = c_estimated(2)/c_estimated(3);
    t_heater_actual = c_actual(2)/c_actual(3);
    c_error_heat(i,1) = t_heater_estimated/t_heater_actual*100;
    c_error_heat(i,2) = c_estimated(1)/c_actual(1)*100;
    c_error_heat(i,3) = c_estimated(2)/c_actual(2)*100;
    c_error_heat(i,4) = c_estimated(3)/c_actual(3)*100;
end

figure
plot(50*heater_sf*9/5+32,c_error_heat)
title('Prediction Accuracy for Varying Heater Air Temperatures')
xlabel('Heater Air Temperature [Deg. F]')
ylabel('Parameter Estimation Accuracy [%]')
axis tight; grid on;
legend('T_h_e_a_t_e_r','C_1 (1/(M*c*Req)','C_2 (Mdot / M * T_h_e_a_t_e_r)','C_3 (Mdot / M)')

req_sf = 0.5:0.05:2;

% Heater temperature sensitivity
c_error_req = zeros(length(req_sf),4);

for i = 1:length(req_sf)
    [t, heater_status, room_temp, c_actual, c_estimated] = thermSim(time, setpoint, outdoor_temp, 'On',0, 0, req_sf(i), 1);
    t_heater_estimated = c_estimated(2)/c_estimated(3);
    t_heater_actual = c_actual(2)/c_actual(3);
    c_error_req(i,1) = t_heater_estimated/t_heater_actual*100;
    c_error_req(i,2) = c_estimated(1)/c_actual(1)*100;
    c_error_req(i,3) = c_estimated(2)/c_actual(2)*100;
    c_error_req(i,4) = c_estimated(3)/c_actual(3)*100;

end

figure
plot(req_sf,c_error_req)
title('Prediction Accuracy for Varying House Heat Loss')
xlabel('R_e_q / R_e_q_ _b_a_s_e_l_i_n_e')
ylabel('Parameter Estimation Accuracy [%]')
axis tight; grid on;
legend('T_h_e_a_t_e_r','C_1 (1/(M*c*Req)','C_2 (Mdot / M * T_h_e_a_t_e_r)','C_3 (Mdot / M)')