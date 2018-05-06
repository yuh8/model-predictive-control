function [V_2,I_2] = abc2dq(V_3,I_3,th)
%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:
%   V_3: three phase voltage
%   I_3: three phase current
%   theta: rotor position
%Output
%   V_2: dq voltage
%   I_2: dg current

V_3 = V_3';
I_3 = I_3';
N = length(th);
I_2 = zeros(2,N);
V_2 = zeros(2,N);

% Compute the three phase angles
th = th(:)';
theta = th;
theta(2,:) = th - 2*pi/3;
theta(3,:) = th + 2*pi/3;

for n = 1:N
K = ParkTransform(theta(:,n));
temp = K*V_3(:,n);
V_2(:,n) = temp(1:2,1); 
temp = K*I_3(:,n);
I_2(:,n) = temp(1:2,1); 
end

function K = ParkTransform(theta)
temp = 1/2;
K = 2/3*[cos(theta(1)), cos(theta(2)), cos(theta(3));...
       -sin(theta(1)), -sin(theta(2)), -sin(theta(3));...
       temp, temp, temp];





