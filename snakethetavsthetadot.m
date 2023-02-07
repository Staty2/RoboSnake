%Code to produce theta and theta dot at different n 

A = 60;
W = (5*pi)/6;
sigma = (2*pi)/3;
n = 1; %change this value for each motor
t = linspace(0, 10, 1000);
theta = A.*sin(W.*t + n.*sigma);



ddt_theta = W.*A.*cos(W.*t +n.*sigma);



tiledlayout(1,2);
nexttile;
plot(t, theta, 'r-')
xlabel('Time')
ylabel('Theta')
nexttile;
plot(t, ddt_theta, 'r-')
xlabel('Time')
ylabel('Theta-dot')





