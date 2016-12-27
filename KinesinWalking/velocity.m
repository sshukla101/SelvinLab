% Written by Saurabh Shukla, Selvin lab, UIUC. 
% Last modified Dec 27, 2016.
 
 
 %% Finding out instantaneous velocity. Select thr transformed folder. Change the exposure time accordingly.
    clear all;   
    exposureTime=0.001; % Input the exposure time.
    CodePath=pwd;
    DataPath=uigetdir;
    cd(DataPath);
    FileIn=dir('*.txt');
    
    %% Code will ask for the time coordinates for calculating the velocity. Time co-ordinates need to be 
    % decided beforehand by looking at the trace.
    
    
    k=1;
    for i=1:length(FileIn)
        fid=fopen(FileIn(i).name);
        Input = textscan(fid,'%f%f','CommentStyle','##');
        yInput = Input{1};
        position=yInput(1:2:end);
        x=yInput(2:2:end);
        plot(x,position);
        
        %Number of velocities calculated per plot
        n_vel=input('how many velocities do you want to calculate from this plot? >> ');
        
        for j=1:n_vel
            %Input start and end position for instaneous velocity
            start_end = input('Start position and end position in matrix form i.e. [800 950]>> ');

        
            v(k)=abs((position(start_end(1)) - position(start_end(2)))/ ((start_end(1)-start_end(2))*exposureTime));
            k=k+1;
        end 
    end
    
    v_avg=mean(v);
    v_max=max(v);
    v_min=min(v);
    
    dlmwrite('inst_velocity.txt',v,'delimiter','\t','precision',3)
    
   
    
    prompt = 'What is the sorbitol concentration?';
    sorbitol = input(prompt,'s');
    
    figure(2)
    bar(v)
    hold on 
    plot([0 length(v)],[v_avg v_avg],'r')
    
    
    title(['Instantaneous velocity of kinesin:',sorbitol, ' M sorbitol'])
    xlabel('samples')
    ylabel('velocity (nm/sec)')
    
    
    
    
    
   
%% Estimating moving average and velocity histogram of the tracce. Select the Transformed folder
    clear all;
    CodePath=pwd;
    DataPath=uigetdir;
    cd(DataPath);
    FileIn=dir('*.txt');
    %% Moving average plots; input the exposure time.
    
    %Moving average
    exp_time=0.00788;
    mov_avg_window=10;
    velocities=[];
    velocity_range=0.1;  %Calculate velocity every 100 ms
    frames=ceil(velocity_range/exp_time);
 
    for i=1:length(FileIn)
        f=dlmread(FileIn(i).name);
        x=f(:,1);
        t=(1:1:length(x))*exp_time;
        moving_avg=tsmovavg(x,'s',mov_avg_window,1); % Calculating the moving average.
        
        plot(t,x,'r');
        hold on;
        plot(t,moving_avg,'b','LineWidth',2);
        xlabel('time (s)');
        ylabel('x (nm)');
        saveas(gcf,['moving_avg_',num2str(i),'.png']);  %Saving the trajectory along with mving average.
        clf;
        
        %Velocity histogram
        iterations=floor(length(x)/frames);
        
        start=1;
        stop=frames;
        for j=1:iterations
            inst_vel=  abs((moving_avg(start)-moving_avg(end))/(frames*exp_time));
            velocities=[velocities;inst_vel];
            start=start+frames;
            stop=stop+frames;
        end    
        
    end
    
    prompt = 'What is the sorbitol concentration?';
    sorbitol = input(prompt,'s');
    
    
    figure(1)
    histogram(velocities(velocities<4000),50);
    title(['velocity histogram:',sorbitol,'M sorbitol']);
    xlabel(' velocity (nm/sec)');
    ylabel('frequency');
    
    
   
