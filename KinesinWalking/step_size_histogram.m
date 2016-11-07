%Created by Saurabh Shukla 
%Nov. 6, 2016
%Grainger Library, UIUC

%%Search txt file in the folder
CodePath=pwd;
DataPath=uigetdir;
cd(DataPath);

%%Search .mat files in the folder
FileIn=dir('*.mat');
FileInput=cell(length(FileIn),1);
FileInputName=cell(length(FileIn),1);

for i=1:length(FileIn)
  FileInput{i}=FileIn(i).name;
  FileInputName{i}=strrep(FileIn(i).name,FileType,'');
end

%Extracting the step_size
CompiledDataPos = [];
CompiledDataNeg = [];
for i=1:length(FileIn)
  Step = load(FileInput{i}); Step = Step.Step;
  dataPos = Step.StepSizeStats(Step.StepSizeStats>0)';
  CompiledDataPos = [CompiledDataPos;dataPos];
  dataNeg = Step.StepSizeStats(Step.StepSizeStats<0)';
  CompiledDataNeg = [CompiledDataNeg;dataNeg];
end


