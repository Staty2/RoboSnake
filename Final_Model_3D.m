close all;
clc;
clear;

x0 = 0;
y0 = 0;
z0 = 0;

z_velocity=2;

incline = 30; %degrees
slope_height = 0;

n = 10; % number of links in centimeters (cm)
l = 10; % length of 1 segment (cm)
mu_n = 0.11*cosd(incline); %normal coefficient of friction
mu_t = 0.2*cosd(incline); %tangent coefficient of friction
T = 4;
M = 1; %total mass
m = M/n; %mass of individual link
g = 10;

r = zeros(2, n + 1); % position vectors
r(:, 1) = [x0; y0]; % initial position


r = zeros(3, n + 1); % position vectors
r(:, 1) = [x0; y0; z0]; % initial position

% making the graph
figure();

%lotting the plane
x = [0, 300, 300, 0];
y = [-150, -150, 150, 150];
z = tand(incline) * x;
patch(x, y, z, [0.4660 0.6740 0.1880]);

hold on;
grid on;
axis([0, 300, -150, 150, 0, 300]);
xlabel('x')
ylabel('y')
zlabel('z')

dT = 0.1;
th = zeros(1, n); % angles for each segment
a = n;
x = [0;0;0];


alpha = 0.2;

for i = 1:length(th) % for each segment in the snake
    [th(i), ~] = gait(i, n); % using the gait function to find angles, ignore the second output
end

for t = 0:dT:4
    for i = 1:length(r)
        pos = x;
        velocity = [0; 0; 0];
        propulsive = [0; 0; 0];

        for j = 2:i-1
            pos = pos + l * [cos(sum(th(1:j))); sin(sum(th(1:j))); tand(incline)];
            [~, th_dot] = gait(j, n); % get the derivative of the gait for each segment
            velocity = velocity + l * [(cos(sum(th(1:j))) * th_dot); (sin(sum(th(1:j))) * th_dot); 0];
            
            % calculate propulsive force
            propulsive = propulsive + alpha * l * th_dot * [cos(sum(th(1:j))); sin(sum(th(1:j)));0];

        end

        % calculate unit tangent and normal vectors
        if i < length(r)
            tangent = (r(:, i+1) - r(:, i)) / l;
            normal = [0, -1, 0; 1, 0, 0; 0, 0, 1] * tangent;
        else
            tangent = [0; 0; 0];
            normal = [0; 0; 0];
        end
        
        % calculate tangential and normal components of velocity
        vel_tangent = dot(velocity, tangent);
        vel_normal = dot(velocity, normal);
       
        % calculate tangential and normal components of propulsive force
        prop_tangent = dot(propulsive, tangent);
        prop_normal = dot(propulsive, normal);
        
        % calculate friction forces per time step
        friction_normal = vel_normal * (mu_n * m * g);
        friction_tangent = vel_tangent * (mu_t * m * g);

        % calculate net force and acceleration
        F_net = propulsive - friction_tangent * tangent - friction_normal * normal - (M * g * sind(incline));
        % calculate the velocity of the center of mass (COM)
        a_COM = (friction_tangent * tangent + friction_normal * normal) / M;
        velocity_COM = a_COM .* l .* [cos(cumsum(th)); sin(cumsum(th)); ones(1,n)] .* dT;

        % update position with velocity and friction force
        r(:, i) = pos(:,1);
        x = x + velocity_COM;
    end

    %this plots the snake
    %%%%%%not sure not plotting +not getting 3d axes. Pretty sure i've checked all the documentation for all the graph stuff and should be fine%%%%%%
        
    % Plot the snake
        snake = zeros(n, 3);
    for i = 1:n-1
        snake(i, 1) = plot3([r(1, i), r(1, i+1)], [r(2, i), r(2, i+1)], [r(3, i), r(3, i+1)],'r');
    end
    for i = ceil(n/4):n
        snake(i, 2) = plot3([r(1, i), r(1, i+1)], [r(2, i), r(2, i+1)], [r(3, i), r(3, i+1)],'r');
    end
    view(3); % set the viewing angle to 3D

    title('Time (seconds):',t)

    pause(dT);
    delete(snake);
    r(:, 1:n) = r(:, 2:n+1);

    a = mod(a+1, n);
    th(1) = th(1) + th(2);
    th(2:n-1) = th(3:n);
    [th(n), ~] = gait(a, n); % update the last angle using the gait function
end

%%%%%%V=(pos(2, end)-y0)/T

function [th, th_dot] = gait(i, n)
    A = pi / 6; % maximum angle for a segment
    W = pi/5; % temporal freq (oscillations per sec)
    t0 = i;
    sigma = pi/3;

    th = A * sin(W * t0 + n * sigma);
    th_dot = A * W * cos(W * t0 + n * sigma);
end

