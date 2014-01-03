function [c_actual, c_estimated, 


%
% Name: househeat_data.m 
% 
% Description:
%
%       Run simulation of house heater and thermostat including a
%       learning mode that identifies the house and heater thermal properties
%       after applying heat for a period of time and observing the response.
%
% Greg Fiore 2/10/2012

%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%

time = 1:1/6:40;                                    % Time steps for simulation [hrs]
setpoint = [time', 70*ones(length(time),1)];        % Set point for normal operations [deg F]
outdoor_temp = [time', 15*sin(2*pi/24*time') + 60]; % Outdoor temperature [deg F]

print_debug = 1;  % Prints data to screen during simulation
plot_debug = 1;    % Plots data for learning algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Learning mode parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

learning_mode = 'Off';   % Status of the learning algorithm ('On' or 'Off')
on_delay = 0.5;         % Time from start of simulation to when heat is turned on [hr]
off_delay = 1;          % Time from start of simulation to when heat is turned off [hr]
temp_threshold = 1;     % Change in room temperature threshold for when heater effects first occur [dec. C]


%%%%%%%%%%%%%%%%%%
% House Geometry %
%%%%%%%%%%%%%%%%%%

lenHouse = 30;        % House length = 30 m
widHouse = 10;        % House width = 10 m
htHouse = 4;          % House height = 4 m
pitRoof = 40/180/pi;  % Roof pitch = 40 deg
numWindows = 6;       % Number of windows = 6
htWindows = 1;        % Height of windows = 1 m
widWindows = 1;       % Width of windows = 1 m
windowArea = numWindows*htWindows*widWindows;   % Total area of windows
wallArea = 2*lenHouse*htHouse + 2*widHouse*htHouse + ...
           2*(1/cos(pitRoof/2))*widHouse*lenHouse + ...
           tan(pitRoof)*widHouse - windowArea;  % Total wall (non-window) area

%%%%%%%%%%%%%%%%%
% Thermal model %
%%%%%%%%%%%%%%%%%

% Heat Input = Mdot * c * (Heater Air Temp - House Temp)
% Room Temp = Integral ( 1/(M*c) * (Heat Input - 1/Req * (Room Temp - Outside Temp))

kWall = 0.038*3600;                     % k is in units of J/sec/m/C - convert to J/hr/m/C multiplying by 3600
LWall = .2;                             % Insulation material in the walls, 0.2 m thick
RWall = LWall/(kWall*wallArea);         % Thermal resistance of walls

kWindow = 0.78*3600;                    % k is in units of J/sec/m/C - convert to J/hr/m/C multiplying by 3600
LWindow = .01;                          % Glass windows, 0.01 m thick
RWindow = LWindow/(kWindow*windowArea); % Thermal resistance of windows

Req = RWall*RWindow/(RWall + RWindow);% Equivalent thermal resistance of the whole house

c = 1005.4;                             % c = cp of air (273 K) = 1005.4 J/kg-K

THeater = 50;                           % Constant temperature of heater air (deg. C)
Mdot = 3600;                            % Air flow rate Mdot = 1 kg/sec = 3600 kg/hr

densAir = 1.2250;                       % Density of air at sea level = 1.2250 kg/m^3
M = (lenHouse*widHouse*htHouse+tan(pitRoof)*widHouse*lenHouse)*densAir;     % Total mass of air in house

TinIC = (70-32)*5/9;                    % Initial indoor temperature (dec. C)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Learning mode parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c_estimated = zeros(3,1); % Initialize estimated parameters 
% c(1):  Estimation of heat loss term c1 = 1/(M*c*Req)
% c(2):  Estimation of heater input c2 = Mdot/M * T_Heater
% c(3):  Estimation of heat transfer c3 = Mdot/M

% Actual parameters (used only for algorithm accuracy compariosns
c_actual = [1/(M*c*Req);...
            Mdot/M * THeater;
            Mdot/M];   

therm_model2;  % Open the model

% Set the switch based on the learning mode
switch learning_mode
    case 'On'
        set_param('therm_model2/Thermostat/Manual Switch','sw','0');
        % Utilize learning mode to turn on the heat and perform system
        % identification algorithm
    otherwise
        set_param('therm_model2/Thermostat/Manual Switch','sw','1');
        % Keep the thermostat functioning normally
end

if print_debug
    fprintf('\n\n\n\nSimulation:  therm_model2\n')
    fprintf('---------------------------------\n')
    fprintf('Learning Mode: %s\n',learning_mode)
    fprintf('Initial House Temp = %3.1f Deg. F\n',TinIC*9/5+32)
    if strcmp(learning_mode,'Off')
        fprintf('House Setpoint = %3.1f Deg. F\n', setpoint(1))
    end
    fprintf('Outdoor Temp = %3.1f to %3.1f Deg. F\n',min(outdoor_temp(:,2)), max(outdoor_temp(:,2)))
    fprintf('Simulation duration = %3.1f hrs\n',time(end))
    fprintf('\n\nRunning Simulation...\n\n')
end


%%%%%%%%%%%%%%%%%%%%%%
% Run the Simulation %
%%%%%%%%%%%%%%%%%%%%%%

sim('therm_model2', [0 time(end)]);

if print_debug
    fprintf('Simulation complete.\n\n')
end

% Extract data returned from the simulation (would be stored by the thermostat in real life)
t = therm_data.time;                                    % Simulation time vector [hr]
heater_status = therm_data.signals.values(:,1);         % Heater status On/Off [1 or 0]
outside_temp = (therm_data.signals.values(:,2)-32)*5/9; % Recorded outside temperature [deg C]  (from wifi weather data)
room_temp = therm_data.signals.values(:,3);             % Recorded room temperature over duration of simulation [dec C]

if plot_debug
    
    %%%%%%%%%%%%%%%%%%%%%
    % Plot thermal data %
    %%%%%%%%%%%%%%%%%%%%%
    
    fh = figure;   set(fh,'Position',[468   338   560   420])     
    ax(1) = subplot(211); 
    plot(t, outside_temp*9/5 + 32, t, room_temp*9/5 + 32);  % Plot outside temperature and house temp
    grid on; axis tight; 
    legend('Outdoor Temp','House Temp')
    title({'Simulation Results';['Learning Mode = ',learning_mode]})
    ylabel('Degrees F')
    
    ax(2) = subplot(212);
    plot(t,heater_status)
    set(ax(2),'YLim',[-0.5, 1.5], 'YTick',[0, 1],'YTickLabel',{'Off','On'})
    grid on
    xlabel('Elapsed Time [hr]')
    legend('Heater Status')
    
    linkaxes(ax,'x')
end



switch learning_mode
    case 'On'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Perform the identification algorithm using the  %
        % "therm_data" output from the simulation         %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Find the indices of significant events in the sequence
        max_id = 1;
        first_id = 1;
        last_id = 1;
        
        % Find the index of the maximum temperature  (could use "find" function, likely not available on platform)
        for i = 2:length(t)
            if room_temp(i) >= room_temp(max_id)
                max_id = i;
            end
        end
        
        % Find the index of the beginning of the temperature rise
        for i = 2:length(t)
            if (room_temp(i) - room_temp(i-1)) >= temp_threshold
                first_id = i;
                break
            end
        end
            
        % Find the index of the end of the cool-down
        t_last = t(max_id) + (off_delay - on_delay);
        for i = 2:length(t)
            if t(i) > t_last
                last_id = i;
                break
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Measurements of c(1) (Heat Loss Coefficient) %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        c_data = zeros(length(t)-max_id-1,1);
        
        for i = 2:length(c_data)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Compute the estimated value of 1/(M*c*Req) for each data point %
            %                                                                %
            %   dT_room         1                                            %
            %   -------  = ----------- * (T_room - T_outside)                %
            %      dt      M * c * Req                                       %
            %                                                                %
            %                 1          dT_room            1                %
            %   C_data = ----------- = - ------- * --------------------      %
            %            M * c * Req       dt      (T_room - T_outdise)      %
            %                                                                %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            idx = max_id + i;
            c_data(i-1) = - ( room_temp(idx) - room_temp(idx-1) ) / ( t(idx) - t(idx-1) ) / ( room_temp(idx-1) - outside_temp(idx-1) );
        end
        
        % Mean value
        c_estimated(1) = mean(c_data);
        c_actual(1) = 1/(M*c*Req);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Measurements of c(2) and c(3) - Heater parameters %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        x = zeros(max_id - first_id-1,1);  % Measurements of T_Room
        y = x;                             % Measurements of C(2) - C(3) * T_room
        
        for i = first_id:max_id
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                                                                             %
            %  dT_room      M * c                             1                           %
            %  ------- =  -------- (T_Heater - T_room) - ----------- (T_room - T_outside) %
            %     dt      Mdot * c                       M * c * Req                      %
            %                                                                             %
            %  dT_room                                                                    %
            %  ------- = C(2) * T_Heater - C(3) * T_room - C(1) * (T_room - T_outside)    %
            %     dt                                                                      %
            %                                                                             %
            %   dT_room                                                                   %
            %  ------- + C(1) * (T_room - T_outside) = C(2) - C(3) * T_room               %
            %     dt                                                                      %
            %                                                                             %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            y(i-first_id+1) = ( room_temp(i) - room_temp(i-1) ) / ( t(i) - t(i-1) ) + c_estimated(1) * ( room_temp(i) - outside_temp(i) );
            x(i-first_id+1) = room_temp(i);
        end
        
        p = polyfit(x,y,1);         % Perform least-squares fit on the data
        
        if plot_debug
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Plot data used for analysis %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fh1 = figure;   set(fh1,'Position',[1032         439         565         371])
            plot(t,room_temp)
            axis tight;
            ylim = get(gca,'YLim');
            hold on
            plot(t(first_id)*[1, 1],ylim,'g--');    % plot start of window
            plot(t(max_id)*[1, 1],ylim,'k--');      % plot maximum temp
            plot(t(last_id)*[1, 1],ylim,'r--');     % plot end of window
            set(gca,'XLim',[on_delay - 0.5*(off_delay-on_delay), off_delay + 0.5*(off_delay-on_delay)])
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Plot thermal coefficient fit %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            fh2 = figure;   set(fh2,'Position',[1032          35         563         327])
            plot(x,y,'.',x,polyval(p,x),'r')
            title('Heat Transfer Coefficient Data')
            grid on
            axis tight;
            xlabel('Room Temperature')
            ylabel('Coefficients')
        end
        
        
        c_estimated(3) = -p(1);    %Mdot/M
        c_estimated(2) = p(2);    % Mdot/M * THeater
        
        c_actual(3) = Mdot/M;
        c_actual(2) = Mdot/M * THeater;
        
        THeater_estimated = c_estimated(2)/c_estimated(3);
        
        if print_debug
            fprintf('---------------------------------------------------\n')
            
            fprintf('Parameter           Actual  Measured  Error [%%]\n')
            fprintf('---------------------------------------------------\n')
            for i = 1:3
                fprintf('C%d                 %6.2f   %6.2f   %6.2f   \n\n',i,c_actual(i), c_estimated(i), (c_estimated(i)-c_actual(i))/c_actual(i)*100)
            end
            fprintf('T_Heater [deg C]   %6.2f   %6.2f   %6.2f   \n\n',THeater, THeater_estimated, (THeater_estimated-THeater)/THeater * 100)
            fprintf('---------------------------------------------------\n')
            
        end
end