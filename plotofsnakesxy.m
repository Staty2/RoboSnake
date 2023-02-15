%set initial y and x co-ord
xin = 1;
yin = 1;
n = 10;
x = zeros(n);
y = zeros(n);
x(1)=xin;
y(1)=yin;

A = (60*pi)/180;%Amplitude
W = (5*pi)/6;%Temporal freq
sigma = (2*pi)/3;%spacial freq
tinitial = 2; %initial time step

d= 0.5;
tval = linspace(1,100, 1000);


for t = 1:length(tval);
    for k = 2:n;
        theta = A.*sin( W.*tval(t) + (k-1).*sigma);
        x(k) = x(k-1) + d.*(sin(theta));
        y(k) = y(k-1) + d.*(cos(theta));

    end
    

plot (y, x, 'r-')
xlim([0 10])
ylim([0 5])
grid on 
pause(0.1)

end
