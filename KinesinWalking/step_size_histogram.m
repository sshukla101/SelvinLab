    %%Search txt file in the folder
    CodePath=pwd;
    DataPath=uigetdir;
    cd(DataPath);
    FileIn=dir('*.mat');    
    %% Extracting the step_size
    CompiledDataPos = [];
    CompiledDataNeg = [];
    for i=1:length(FileIn)
        Step = load(FileIn(i).name); Step = Step.Step;
        dataPos = Step.StepSizeStats(Step.StepSizeStats>0)';
        CompiledDataPos = [CompiledDataPos;dataPos];
        dataNeg = Step.StepSizeStats(Step.StepSizeStats<0)';
        CompiledDataNeg = [CompiledDataNeg;dataNeg];
    end
    %% Making the histograms of the data
    prompt = 'What is the sorbitol concentration?';
    sorbitol = input(prompt,'s');
    
    figure(1);
    histogram([CompiledDataPos; CompiledDataNeg],100);
    xlim([-40 40]);
    title(['Step size histogram:',sorbitol, ' M sorbitol'])
    xlabel('step size (nm)')
    ylabel('frequency')

    
%% Finding out instantaneous velocity
    exposureTime=0.00788;
    CodePath=pwd;
    DataPath=uigetdir;
    cd(DataPath);
    FileIn=dir('*.txt');
    for i=1:length(FileIn)
        fid=fopen(FileIn(i).name);
        Input = textscan(fid,'%f%f','CommentStyle','##');
        yInput = Input{1};
        position=yInput(1:2:end);
        x=yInput(2:2:end);
        plot(x,position);
        
        %Input start and end position for instaneous velocity
        prompt = 'Start position';
        start_l = input(prompt);
        prompt = 'End position';
        end_l = input(prompt);
        
        v(i)=abs((position(start_l) - position(end_l))/ ((start_l-end_l)*exposureTime));
    end
    
    v_avg=mean(v);
    v_max=max(v);
    v_min=min(v);
    
     
        
        
        
        
        
        
            
    
