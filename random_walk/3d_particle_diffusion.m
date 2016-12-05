%MATLAB code for 3d diffusion of particles.
%File creates a output file where columns represent [x y z timestep particle_number]
%It can generate the animation of the trajectories of particles
%Author: Saurabh Shukla, Dec 5, 2016. 


clear all;

N=100;          %number of particles
Nsteps=10000;    %Number of steps
tau= 0.1;    %Time step 
time=100;    %total simulation time.
D=1;       % diffusion coefficient

k= sqrt(2*D*tau);  %scaling factor


x_0= 10*rand(N,1);
y_0= 10*rand(N,1);
z_0= 10*rand(N,1);


x=  k * randn(N,Nsteps);  %each row of x is journey of a particle at different times.
y=  k * randn(N,Nsteps);
z=  k * randn(N,Nsteps);
    

for i=1:N
    x(i,:)=x(i,:)+x_0(i);
    y(i,:)=y(i,:)+y_0(i);
    z(i,:)=z(i,:)+z_0(i);
end

%% Saving the data
traj=zeros(N*Nsteps,5);
time=tau*linspace(1,Nsteps,Nsteps);


output=zeros(N*Nsteps,5);

j=1;
for i=1:N
    output(j:i*Nsteps,1)=x(i,:)';
    output(j:i*Nsteps,2)=y(i,:)';
    output(j:i*Nsteps,3)=z(i,:)';
    output(j:i*Nsteps,4)=time';
    output(j:i*Nsteps,5)=i*ones(Nsteps,1);
    j=j+Nsteps;
end

save 3d_diffusion.mat output


%% Trajectory simulation of data
for i=1:Nsteps
    c = linspace(1,10,N);
    plot3(x(:,i),y(:,i),z(:,i))
    axis([0 10 0 10 0 10]);
    hold on
    pause(0.1);
end


%% Scatter simulaton of the data


for i=1:Nsteps
    c = linspace(1,10,N);
    scatter3(x(:,i),y(:,i),z(:,i),30,c,'filled')
    axis([0 10 0 10 0 10]);
    drawnow
    pause(0.1);
end
    

    

    

