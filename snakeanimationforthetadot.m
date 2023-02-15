%Calculates theta dot between motors
%All values same as the other code
A = 60;
W = (5*pi)/6;
sigma = (2*pi)/3;

tinitial = 2;
nval = linspace (0,100, 100);

theta_atn = W.*A.*cos( W.*tinitial + nval.*sigma);



plot(nval, theta_atn,  'r-')
tval = linspace(1,100,100);

for t= 1:length(tval);
    theta = W.*A.*cos( W.*tval(t) + nval.*sigma);
    plot(nval,theta,'o')
    ylabel('theta-dot')
    xlabel('motor number')

    pause(.1)

end
