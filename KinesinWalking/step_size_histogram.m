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

    
