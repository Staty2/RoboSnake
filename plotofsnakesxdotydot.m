
%set initial y and x co-ord
xin = 1;
yin = 1;
n = 10;
xdot = zeros(n);
ydot = zeros(n);
xdot(1)=2;
ydot(1)=2;

x= zeros(n);
y = zeros(n);
x(1)=xin;
y(1)=yin;

xpos = zeros(n);
ypos = zeros(n);

A = (60*pi)/180;%Amplitude
W = (5*pi)/6;%Temporal freq
sigma = (2*pi)/3;%spacial freq


d= 0.5;
tval = linspace(1,100, 100000);
for k = 2:n;
    theta = A.*sin( W.*tval(1) + (k-1).*sigma);
    thetadot = W.*A.*cos( W.*tval(1) + (k-1).*sigma);
    
    xdot(k) = xdot(k-1) + d.*(sin(theta)).*thetadot;
    ydot(k) = ydot(k-1) + d.*(cos(theta)).*thetadot;
    x(k) = x(k-1) + d.*(sin(theta));
    y(k) = y(k-1) + d.*(cos(theta));
    
end
sumy = 0;
sumx = 0;
for s = 1:n;
    sumx = sumx + xdot(s);
    sumy = sumy + ydot(s);
end;
xdot(1) = sumx;
ydot(1) = sumy;


for t = 2:length(tval);
    for k = 1:n;
        xpos(k) = x(k) + xdot(k).*0.001;
        ypos(k) = y(k) + ydot(k).*0.001;
    end




    for k = 2:n;
        theta = A.*sin( W.*tval(t) + (k-1).*sigma);
        thetadot = W.*A.*cos( W.*tval(t) + (k-1).*sigma);
        
        xdot(k) = xdot(k-1) + d.*(sin(theta)).*thetadot;
        ydot(k) = ydot(k-1) + d.*(cos(theta)).*thetadot;
        x(k) = xpos(k-1) + d.*(sin(theta));
        y(k) = ypos(k-1) + d.*(cos(theta));
        
    end;
    sumy = 0;
    sumx = 0;
    for s = 1:n;
        sumx = sumx + xdot(s);
        sumy = sumy + ydot(s);
    end;
    xdot(1) = sumx;
    ydot(1) = sumy;


    plot (ypos, xpos, 'r-')
    xlim([-30 30]);
    ylim([-30 30]);
    pause(0.1)
end
