%set initial y and x co-ord
xin = 1;
yin = 1;
n = 100;
x = zeros(n);
y = zeros(n);
x(1)=xin;
y(1)=yin;
xdot = zeros(n);
ydot = zeros(n);
xdot(1)=0;
ydot(1)=1;


A = (60*pi)/180; %Amplitude
W = (5*pi)/6;%Temporal freq
sigma = (2*pi)/9;%spacial freq


d= 0.2;
tval = linspace(1,100, 10000);


ntrotation = [0;0]
vn = zeros([n-1, length(tval)]);
vt = zeros([n-1, length(tval)]);
vnvt = [0;0]
for k = 2:n
        theta = A.*sin( W.*tval(1) + (k-1).*sigma);
        x(k) = x(k-1) + d.*(sin(theta));
        y(k) = y(k-1) + d.*(cos(theta));

end

for t = 1:length(tval)
    %calculating positions for x and y based of previous equations
    
    for k = 2:n
        theta = A.*sin( W.*tval(t) + (k-1).*sigma);
        thetadot = W.*A.*cos(W.*tval(t) + (k-1).*sigma);
        xdot(k) = xdot(k-1) + d.*(cos(theta)).*thetadot;
        ydot(k) = ydot(k-1) - d.*(sin(theta)).*thetadot;
       



    end

    

   for k = 2:n
       x(k) = x(k) + xdot(k).*0.01;
       y(k) = y(k) + ydot(k).*0.01;
   end
   
   



plot (y([2:n]), x([2:n]), 'r-')


ylabel('snake moving in y direction')
xlabel('snake moving in x direction')
xlim([0 10])
ylim([0 5])
title('Time:',tval(t) )
grid on
pause(0.01)
end
