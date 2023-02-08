A = 60;%Amplitude
W = (5*pi)/6;%Temporal freq
sigma = (2*pi)/2;%spacial freq
%Interesting when you divide the frequencies by different numbers
tinitial = 2; %initial time step
%number of motors can be changed to generate the point masses
nval = linspace (0,1000, 1000);

theta_atn = A.*sin( W.*tinitial + nval.*sigma);


tiledlayout(1,1)
plot(nval, theta_atn,  'r-')
tval = linspace(1,100,10000); %Change time values and spaces between them 

for t= 1:length(tval);
    theta = A.*sin( W.*tval(t) + nval.*sigma);%equation
    plot(nval,theta,'r-')
    xlabel('motor number')
    ylabel('theta')
    pause(.1)
end




