close all;
clc;
clear;

%%%%i don't know how to turn this into 3d plot, i was using plot3%%%%%%

x0 = 0;
y0 = 0;
z0 = 0;

z_velocity=10;

n = 10;
l = 10; % length of 1 segment
mu_n = 0.11; %  normal kinetic friction coefficient
mu_t = 0.2; %  tangential kinetic friction coefficient
alpha = (5*pi)/180;
T=4;

r = zeros(3, n + 1); % position vectors
r(:, 1) = [x0; y0; z0]; % initial position

% making the graph
figure();
hold on;
grid on;
axis([0, 200, -100, 100, 0, 10]);

ylabel('snake moving in y direction')
xlabel('snake moving in x direction')
zlabel('snake moving in z direction')

dT = 0.1;
th = zeros(1, n); % angles for each segment
a = n;
x = [0;0;0];


for i = 1:length(th) % for each segment in the snake
    [th(i), ~] = gait(i, n); % using the gait function to find angles, ignore the second output
end

for t = 0:dT:T
    for i = 1:length(r)
        pos = x;
        velocity = [0; 0; z_velocity];

        for j = 2:i-1
            pos = pos + l * [cos(sum(th(1:j))); sin(sum(th(1:j))); sin(alpha)*cos(sum(th(1:j)))];
            [~, th_dot] = gait(j, n); % get the derivative of the gait for each segment
            velocity = velocity + l * [(cos(sum(th(1:j))) * th_dot); (sin(sum(th(1:j))) * th_dot); z_velocity/l];

        end

        % calculate unit tangent and normal vectors 
        if i < length(r)
            tangent = (r(:, i+1) - r(:, i)) / l;
            %%%%%%%don't understand this normal bit well enough to know
            %%%%%%%what to times the z part by
            normal = [0, -1, 0; 1, 0, 0; 1, 0, 0] * tangent;
        else
            tangent = [0; 0; 0];
            normal = [0; 0; 0];
        end

        % calculate tangential and normal components of velocity
        vel_tangent = dot(velocity, tangent);
        vel_normal = dot(velocity, normal);
        
        %friction forces per time step
        friction_normal = vel_normal.*(mu_n*n*l);
        friction_tangent = vel_tangent.*(mu_t*n*l);
        gravity_force = m %I have hit the point where my brain is confused by gravity so I am going to bed
      
        %%not sure we need this section??%%
        velocity_COM = sum(velocity(1,:))/n;

        % update position with velocity and friction force
        r(:, i) = pos;
        x = x + [velocity_COM;0]*dT;
    end

    %this plots the snake
    %%%%%%not sure not plotting +not getting 3d axes. Pretty sure i've checked all the documentation for all the graph stuff and should be fine%%%%%%
    snake = zeros(n, 3);
    for i = 1:n-1
        snake(i, 1) = line([r(1, i), r(1, i+1)], [r(2, i), r(2, i+1)], [r(3, i), r(3, i+1)]);
    end
    for i = ceil(n/4):n
        snake(i, 2) = line([r(1, i), r(1, i+1)], [r(2, i), r(2, i+1)], [r(3, i), r(3, i+1)]);
    end

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

