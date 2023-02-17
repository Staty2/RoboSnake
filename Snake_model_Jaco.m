close all;
clc;
clear;

x0 = 0;
y0 = 0;

n = 10;
l = 10; % length of 1 segment
mu_n = 0.11;% kinetic friction coefficient
mu_t = 0.2;

r = zeros(2, n + 1); % position vectors
r(:, 1) = [x0; y0]; % initial position

% making the graph
figure();
hold on;
grid on;
axis([0, 200, -100, 100]);
ylabel('snake moving in y direction')
xlabel('snake moving in x direction')

dT = 0.1;
th = zeros(1, n); % angles for each segment
a = n;
x = [0;0];

alpha = 0.2;
beta = 0.4;

for i = 1:length(th) % for each segment in the snake
    [th(i), ~] = gait(i, n); % using the gait function to find angles, ignore the second output
end

for t = 0:dT:4
    for i = 1:length(r)
        pos = x;
        velocity = [0; 0];

        for j = 2:i-1
            pos = pos + l * [cos(sum(th(1:j))); sin(sum(th(1:j)))];
            [~, th_dot] = gait(j, n); % get the derivative of the gait for each segment
            velocity = velocity + l * [(cos(sum(th(1:j))) * th_dot); (sin(sum(th(1:j))) * th_dot)];
            
        end

        % calculate unit tangent and normal vectors
        if i < length(r)
            tangent = (r(:, i+1) - r(:, i)) / l;
            normal = [0, -1; 1, 0] * tangent;
        else
            tangent = [0; 0];
            normal = [0; 0];
        end
        
        % calculate tangential and normal components of velocity
        vel_tangent = dot(velocity, tangent);
        vel_normal = dot(velocity, normal);
        
        %friction forces per time step
        friction_normal = vel_normal.*(mu_n*n*l);
        friction_tangent = vel_tangent.*(mu_t*n*l);
     
        velocity_COM = sum(velocity(1,:))/n;

        % update position with velocity and friction force
        r(:, i) = pos;
        x = x + [velocity_COM;0]*dT;
    end
    
    %this plots the snake
    snake = zeros(n, 2);
    for i = 1:n-1
        snake(i, 1) = line([r(1, i), r(1, i+1)], [r(2, i), r(2, i+1)]);
    end
    for i = ceil(n/4):n
        snake(i, 2) = line([r(1, i), r(1, i+1)], [r(2, i), r(2, i+1)]);
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

function [th, th_dot] = gait(i, n)
    A = pi / 6; % maximum angle for a segment
    W = pi/5; % temporal freq (oscillations per sec)
    t0 = i;
    sigma = pi/3;
   
    th = A * sin(W * t0 + n * sigma);
    th_dot = A * W * cos(W * t0 + n * sigma);
end

