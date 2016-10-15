function [ angle_change ] = delta( ang_vel, ang_accel, t)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
deltaX = ang_vel(1,1) * t + ang_accel(1,1)*t^2;
deltaY = ang_vel(1,2) * t + ang_accel(1,2)*t^2;
deltaZ = ang_vel(1,3) * t + ang_accel(1,3)*t^2;
angle_change = [ deltaX deltaY deltaZ];
end

