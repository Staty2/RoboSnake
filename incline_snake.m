%variables
theta_0 =(30*pi)/180; %start direction in RADIANS
%for some angles it goes backwards but not sure why??
N= 11; %number of links
r_00=[1,2,0]'; %position in x,y,z currently not on any slope
d=0.1; %10cm links
T=5; %max time
a=(5*pi)/180; %slope angle
M=10; %mass of snake
m=10/N; %mass of each node (assuming constant mass)

theta_list = zeros(1,N);
for i=1:2:N
    theta_list(i)=theta_0;
    theta_list(i+1)=-theta_0;
end

%intial pos of first node (end of snake)
initial_position = zeros(3,N);
initial_position(:,1) = r_00;

%intial pos of whole snake
for i=1:N-1
   initial_position(1,i+1)=initial_position(1,i)+d*sin(theta_list(i));
   initial_position(2,i+1)=initial_position(2,i)+d*cos(theta_list(i));
   initial_position(3,i+1)=initial_position(3,i)+sin(a)*d*cos(theta_list(i));
  
end

%this is shit naming. The 0.1 is how you set the speed - a sort of
%frequency
time_span=linspace(1,T,T/0.1);
p=length(time_span);

position_matrix=zeros(3,N,p);
position_matrix(:,:,1)=initial_position;
position=initial_position;



mass_position=zeros(3,p);

for j=1:p

    %cofm position calc for each time step
    x_val=position(1,:);
    y_val=position(2,:);
    z_val=position(3,:);
  
    x_mass_i=m.*x_val;
    y_mass_i=m.*y_val;
    z_mass_i=m.*z_val;

    mass_position(1,j)=(1/M)*sum(x_mass_i);
    mass_position(2,j)=(1/M)*sum(y_mass_i);
    mass_position(3,j)=(1/M)*sum(z_mass_i);

    %plotting x,y positions for each time step
    plot3(x_val,y_val,z_val)
    xlim([0,2])
    ylim([0,10])
    zlim([0,0.5])
    pause(0.1)

    %changing start noded position adn thetas for new snake position
    position(:,1) = position(:,2);
    theta_list=-theta_list;


    for i=1:N-1
       position(1,i+1)=position(1,i)+d*sin(theta_list(i));
       position(2,i+1)=position(2,i)+d*cos(theta_list(i));
       position(3,i+1)=position(3,i)+sin(a)*d*cos(theta_list(i));
       
      
    end

    position_matrix(:,:,j+1)=position;
end

mass_position;
position_matrix;


%v simple calc assuming constant speed and no forces:
end_pos=position(:,1);
d=norm(end_pos)-norm(r_00);
v=d/T  %not realistic this is a v slow snake


