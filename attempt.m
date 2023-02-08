
x0 = 2;
y0 = 2;
n = 5;
l = 2;
T = linspace(1,100,10000); %Change time values and spaces between them 

r = zeros(2,n); % position vectors;
r(:,1)=[x0;y0]; % initial position

x = [];

position = [x0;y0];

for t = 1:length(T)

    for i = 2:n
   
    position = position + l*[sum(cos(gait(i,T(t))));sum(sin(gait(i,T(t))))];
    
    r(:,i) = position;
    end
    r

    plot(r(1,:),r(2,:))
    xlim([-50 50])
    ylim([-50 50])
    pause(0.5)
    
end




function theta = gait(nval,tval)

A = 60;%Amplitude
W = (5*pi)/6;%Temporal freq
sigma = (2*pi)/3;%spacial freq
%Interesting when you divide the frequencies by different numbers

theta = A.*sin( W.*tval + nval.*sigma);%equation

end

