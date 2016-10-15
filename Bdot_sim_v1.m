clc, clear all
addpath('Orbit3D')
addpath('math')
addpath('igrf')
format bank
disp(['Simulation started at ' datestr(now)])
tic %Initiate runtime counter

%% ***CONSTANTS***
g_parameter = 3.986E14; %graviational parameter of earth, m^3 / s^2
earth_radius = 6378; %radius of earth, km
secs_per_day = 86164; %seconds in a solar day
momInertiaMatrix = [0.0458 0 0; 0 0.0917 0; 0 0 0.1192];
%Invert the mom inertia matrix for use in a loop later
invMomInertiaMatrix = inv(momInertiaMatrix);

%% ***INPUTS***
sim_time = 96; %minutes in simulation, use whole numbers
    %95.63 = 1 orbit
time_step = .05; %seconds between data sample--use numbers below 1
    %DO NOT SET BELOW 0.01 LEST YOU WISH TO MEET THE KRAKEN
% Empty matrices for storing values later
k = -250000; %B-dot gain value for de-tumble
B_field_orb = zeros(sim_time/time_step,3);
B_field_body = zeros(sim_time/time_step,3);
mag_moment = zeros(sim_time/time_step,3); 
torque = zeros(sim_time/time_step,3);
ang_mom = zeros(sim_time/time_step,3);
ang_accel = zeros(sim_time/time_step,3);
error_quat = zeros(sim_time/time_step,3);
ang_vel(1,:) = [8 -6 4];%Initial angular velocity at deployment, dps
meas_quat = [.707 0 0 .707]; %Initial attitude
orb_alt = 550; %mean altitude of orbit, km
init_pos = [0 0]; %initial position, [latitude longitude]
inclin = 0; %orbital inclination, degrees
    %0 < inclin < 89 means the satellite rotates with earth
    %91 < inclin < 179 means the satellite rotates against earth
eccent = 0; %eccentricity of orbit, 0 is circular
RA_ascend = 0; %right ascension of ascending node
arg_per   = 0; %argument of perigee
true_anom = 0;  %true anomaly of departure
semi_major = orb_alt + earth_radius; %semi_major axis height, km


%% ***CALCULATE***
orb_vel = sqrt(g_parameter/(semi_major * 1000)); %orbital velocity, m/s
period = (2*pi*(semi_major * 1000)) / orb_vel;
num_orb = (sim_time * 60) / period;
fprintf('Simulating %2.2f Orbits \n', num_orb)
time = datenum(datetime);
time_of_day = mod(time, secs_per_day); % seconds in current day
gha = (time_of_day)/secs_per_day*360; %greenwich hour angle, degrees
% Find the latitude and longitude vectors from orbit specs
[lat, long] = Orbit3D(RA_ascend, arg_per, true_anom, inclin, semi_major,...
    eccent, time_step, num_orb);
% Call IGRF to create a matrix of magnetic field values in the NED frame
[Bx_NED, By_NED, Bz_NED] = igrf(time, lat, long, orb_alt);
B_field_NED = [Bx_NED, By_NED, Bz_NED];
% Convert magnetic field vector from NED to orbital frame, row by row
% Creates a matrix of mag field components in orbital frame
c = 1; %Counter for next for loop for matrix indices
for i = 0:(time_step/60):sim_time
    B_field_orb(c,:) = ned_to_orb(inclin, period, i*60)*B_field_NED(c,:)';
    c = c+1;
end
% Loop for each B field vector, converting to body frame at the start
for j = 1:size(B_field_orb, 1)
    % Calculate current vectors based off the orbital B field and inital
    % angular velocity and attitude defined in INPUTS
    meas_quat(j,:) = meas_quat(j,:)/norm(meas_quat(j,:));
    B_field_body(j,:) = v_rot_q(B_field_orb(j,:)', meas_quat(j,:)')*10^-9;
    mag_moment(j,:) = k*cross(B_field_body(j,:),ang_vel(j,:));
    torque(j,:) = cross(mag_moment(j,:), B_field_body(j,:));
    ang_mom(j,:) = momInertiaMatrix * ang_vel(j,:)';
    %ang_accel(j,:) = calc_ang_accel(torque(j,:), ang_vel(j,:), ang_mom(j,:), momInertiaMatrix);
    ang_accel(j,:) = (invMomInertiaMatrix*torque(j,:)')';
    % Find theta_i - theta_f; reference the delta function for more info
    angle_change = delta(ang_vel(j,:), ang_accel(j,:), time_step);
    % Convert the change in angle to a quaternion to combine with the 
    % measured quaternion value
    change_quat = mat2quat(dcm(angle_change));
    change_quat = change_quat/norm(change_quat);
    % Update angular velocity, and attitude except for the last loop
    if j < size(B_field_orb, 1)
        meas_quat(j+1,:) = quat_mult2(change_quat, meas_quat(j,:));
        ang_vel(j+1,:) = ang_accel(j,:)*time_step + ang_vel(j,:);
    end
end

%% ***OUTPUT***
% NED frame component graph
total_time = (0:(time_step / 60):sim_time)';
figure()
plot(total_time, B_field_NED)
ylim([-3.5*10^4 3.5*10^4])
title('Magnetic Field (NED frame)')
legend('Bx', 'By', 'Bz')
xlabel('time, minutes')
ylabel('nano-Tesla')
% Orbital frame component graph
figure()
plot(total_time, B_field_orb)
ylim([-3.5*10^4 3.5*10^4])
title('Magnetic Field (orbital frame)')
legend('Bx', 'By', 'Bz')
xlabel('time, minutes')
ylabel('nano-Tesla')
% Angular velocity graph
figure()
plot(total_time, ang_vel)
title('Angular Velocity')
legend('wx',  'wy', 'wz')
xlabel('time, minutes')
ylabel('angular velocity, deg/s')

% Display sizes to make sure they all match; only for testing code
runtime = toc; %end runtime counter
fprintf('Total Simulation Time = %4.2f \n', runtime)