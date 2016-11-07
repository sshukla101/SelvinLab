%% Perform FIONA analysis on all tiff files and save as text file
% Output is a matrix of [PhotonNo xPrecision yPrecision xCenter yCenter RSquare ErrorStatus]
% Still need to work on allowing it to do parallel processing
% Get folder for files and assign filenames to FileInputName

    % Sum frames to lengthen exposure time
    SumFactor = 1;              % Number of frames to be added together

    % List of parameters                              
    ActPix = 16000;             % Actual pixel size in nanometer
    ObjMag = 100;               % Objective magnification
    AddMag = 1.6;               % Additional magnification

    CCDsens = 12;             % CCD sensitivity of the camera at specific readout rate and pre-amp setting. See note below for values
    EMgain = 100;                % Electron multiplying (EM) gain setting of camera during acquisition

    AskPath = 1;                % (0 or 1) Prompt Matlab to start a dialog box asking directory
    DataPath = 'C:\Users\';     % Path used if AskPath is zero, for automation purposes
    FileType = '.tif';          % Input file type
    Parallel = 0;               % (0 or 1) for no and with parallel computing using matlabpool
    
    % Get directories
    CodePath = pwd;
    if AskPath == 1
        DataPath = uigetdir;
    end
    cd(DataPath);

    % Get all the file names and assign them under FileInput
    FileIn = dir(['*' FileType]);
    FileInput = cell(length(FileIn),1);
    FileInputName = cell(length(FileIn),1);
    for ind=1:length(FileIn)
        FileInput{ind}=FileIn(ind).name;
        FileInputName{ind}=strrep(FileIn(ind).name,FileType,'');
    end
    
% Applying FIONA to each file and generate output file
PixelSize = ActPix/(ObjMag*AddMag);     % Pixel size in nanometer
CountToPhoton = CCDsens/EMgain;         % Count to photon conversion
a2 = PixelSize*PixelSize;               % Square of PixelSize

% Allow use of multiple cores if the 'Parallel' option is 1
if Parallel == 1
    % Check if existing matlabpool is opened. Matlabpool allows the use
    % of more than one core of the computer
    IsOpen=matlabpool('size')>0;
    if IsOpen == 0
        matlabpool open;
    end
    % Find out the number of cores in the computer. If the number of cores
    % returned is less than expected, it can be increased by going to
    % Parallel>Manage Configurations, right click on 'local', click 
    % 'Properties', under the 'Scheduler', you can find 'Number of workers
    % available to scheduler (ClusterSize). You can increase this number. You
    % can find out if the computer is making use of all of its cores using task
    % manager in Windows or Activity Monitor in Mac. In my experience, 
    % Matlabpool will only use half of the available cores by default.
    CoreNo = matlabpool('size');
else
    CoreNo = 1;
end

% Looping over all files
for ind=1:length(FileIn)
    % Getting info from image to make it faster to upload multiple frames
    cd(DataPath);
    info = imfinfo(FileInput{ind});
    num_images = numel(info);
    
    FinalFrameNo = floor(num_images/SumFactor);
    %FinalFrame = double(zeros(info(1).Height,info(1).Width,FinalFrameNo));

    % Allocating memory to variables
    %Output = zeros(num_images,7);
    Output = zeros(FinalFrameNo,7);
    
    tic;
    % Loop over all frames
    for k = 0:FinalFrameNo-1
        data = double(zeros(info(1).Height,info(1).Width));
        for j = 1:SumFactor
            data = data+double(imread(FileInput{ind}, (SumFactor*k)+j, 'Info', info));
        end
        %dimension=size(data);               % Find out the dimensions
        [ny,nx]=size(data);                 % Find out the dimensions
        %xAxis = 1:nx;                       % Create the x-Axis gridding
        %yAxis = 1:ny;                       % Create the y-Axis gridding
        grid = [nx ny 1:nx 1:ny];         % Gridding input for gauss2dfunct and gauss2dfit
        tilt=0;
        %tiltVal=0;

        % Change to code directory
        cd(CodePath);

        % Parameters: p(1): z-offset, p(2): amplitude, p(3): xStdev, p(4): yStdev, p(5): xCenter, p(6): yCenter, p(7): tilt.
        try
            [popt,resnorm,residual,ret]=gauss2dfit(data,grid,tilt);

            if popt(5)>0 && popt(5)<nx && popt(6)>0 && popt(6)<ny

                % Getting center and precision
                xCenter = popt(5)*PixelSize;        % Center of x in nanometer
                yCenter = popt(6)*PixelSize;        % Center of y in nanometer

                PixelSize2 = PixelSize*PixelSize;
                sx2 = popt(3)*popt(3)*PixelSize2;   % Square of xStdev
                sx4 = sx2*sx2;                      % xStdev^4
                sy2 = popt(4)*popt(4)*PixelSize2;   % Square of yStdev
                sy4 = sy2*sy2;                      % yStdev^4

                % To estimate b, which is the standard deviation of the background,
                % we'll look at the z-offset (popt(1)) and calculate the standard
                % deviation based on anything below the z-offset. Before that we
                % would want to make everything 3 standard deviations away from
                % our spot to be (z-offset + 1) or NaN.

                % Find the limits of the data
                xmin = floor(popt(5)-4*popt(3));
                xmax = ceil(popt(5)+4*popt(3));
                ymin = floor(popt(6)-4*popt(4));
                ymax = ceil(popt(6)+4*popt(4));

                if xmin < 1; xmin = 1; end
                if xmax > nx; xmax = nx; end
                if ymin < 1; ymin = 1; end
                if ymax > ny; ymax = ny; end

                %zfit = gauss2dfunct(popt,grid);
                %zfit = reshape(zfit,ny,nx);         % Convert 1d array into 2d array

                z = reshape(data,ny*nx,1); 
                SumSquaresTotal = sum(z.*z);
                RSquare = 1-(resnorm/SumSquaresTotal);

                %residual = data - zfit;                % Compute residual
                residual = reshape(residual,ny,nx);     % Extract residual
                residual = -residual;                   % Invert such that residual = data - zfit
                residual(xmin:xmax,ymin:ymax)=0;        % Set values of residual around center to be NaN
                residual = residual(residual<0);
                b = sqrt(sum(residual.*residual)/(length(residual)-1));
                b = b*CountToPhoton;                    % Standard deviation of the background
                b2 = b*b;                               % Square of background
                PhotonNo = abs(2*pi*popt(2)*popt(3)*popt(4)*CountToPhoton);   % Number of Photons calculated using volume under gaussian, which is 2*pi*A*stdev(x)*stdev(y)
                PhotonNo2 = PhotonNo*PhotonNo;          % Square of PhotonNo
                xPrecision = sqrt((sx2/PhotonNo) + (a2/(12*PhotonNo)) + (8*pi*sx4*b2/(a2*PhotonNo2)));
                yPrecision = sqrt((sy2/PhotonNo) + (a2/(12*PhotonNo)) + (8*pi*sy4*b2/(a2*PhotonNo2)));

                ErrorStatus = 0;

                cd(DataPath);
                Output(k+1,:) = [PhotonNo xPrecision yPrecision xCenter yCenter RSquare ErrorStatus];
            else
                cd(DataPath);
                ErrorStatus = 1;
                Output(k+1,:) = [0 0 0 0 0 0 ErrorStatus];
            end
        catch
            cd(DataPath);
            ErrorStatus = 2;
            Output(k+1,:) = [0 0 0 0 0 0 ErrorStatus];
        end
    end
    toc;
    
    if 0
        % Make points with errors equal previous values
        ErrorRow = find(Output(:,7)>0);   % Row at which there are errors
        if isempty(ErrorRow)==0
            % Take out xOut at the beginning
            loop=1;
            while loop == 1
                if ErrorRow(1)==1
                    Output(1,:)=[];
                    ErrorRow(1)=[];
                    ErrorRow=ErrorRow-1;
                else
                    loop=0;
                end
            end
            % Make points on ErrorRow equal previous values
            for l=1:length(ErrorRow)
                Output(ErrorRow(l),1:6)=Output((ErrorRow(l)-1),1:6);
            end
        end
    end

    % Save Output
    dlmwrite([FileInputName{ind} '(' num2str(SumFactor) 'x).txt'],Output,'\t');
end

clear i;
% Notes on CCD sensitivity for different cameras
% Sub-zero (Readout rate on left (e.g. 10 MHz 14 bit), Preamp setting in
% bracket (e.g. 1x, 2.3x, 4.5x), CCD sensitivity after colon)
%   10 MHz 14 bit: 67.42 (1x), 26.85 (2.3x), 12.27 (4.5x)
%    5 MHz 14 bit: 59.86 (1x), 24.36 (2.3x), 10.82 (4.5x)
%    3 MHz 14 bit: 59.84 (1x), 24.17 (2.3x), 10.72 (4.5x)
%    1 MHz 16 bit: 24.33 (1x),  9.74 (2.3x),  4.26 (4.5x)
% 
% Johnny Cage 
% Readout rate on left (e.g. 10 MHz 14 bit), Preamp setting in
% bracket (e.g. 1x, 2.4x, 4.9x), CCD sensitivity after colon
%   10 MHz 14 bit: 64.5 (1x), 26.3 (2.4x), 11.9 (4.9x)
%    5 MHz 14 bit: 55.0 (1x), 23.9 (2.4x), 10.4 (4.9x)
%    3 MHz 14 bit: 55.8 (1x), 23.6 (2.4x), 10.3 (4.9x)
%    1 MHz 16 bit: 23.3 (1x),  9.1 (2.4x),  4.2 (4.9x)
% 
% Scorpion
% Readout rate on left (e.g. 10 MHz 14 bit), Preamp setting in
% bracket (e.g. 1x, 2.3x, 4.5x), CCD sensitivity after colon
%   10 MHz 14 bit: 67.75 (1x), 27.00 (2.5x), 12.13 (5.2x)
%    5 MHz 14 bit: 57.20 (1x), 23.20 (2.5x), 10.33 (5.2x)
%    3 MHz 14 bit: 56.17 (1x), 23.07 (2.5x), 10.05 (5.2x)
%    1 MHz 16 bit: 23.02 (1x),  9.29 (2.5x),  4.18 (5.2x)
%% Plot stepsizes from SICstepfinder analysis (Use t-test files)    
    %% Search txt file in the folder
    CodePath = pwd;
    DataPath = uigetdir;
    cd(DataPath);
    FileType = '.txt';
    FileIn = dir(['*' FileType]);
    FileInput = cell(length(FileIn),1);
    FileInputName = cell(length(FileIn),1);
    for ind=1:length(FileIn)
        FileInput{ind}=FileIn(ind).name;
        FileInputName{ind}=strrep(FileIn(ind).name,FileType,'');
    end
    %% Open, plot file and analyze with SIC Stepfinder
    ind = 1;
    fid=fopen(FileInput{ind});
    Input = textscan(fid,'%f%f','CommentStyle','##');
    yInput = Input{1};
    %plot(yInput);
    %DataPath = pwd;
    %CodePath = 'C:\Users\tjioe2\Documents\MATLAB\Codes\TugOfWar\FIONA';
    cd(CodePath);
    Step=SICstepFinder(yInput');
    [HistStep,xout] = hist(Step.StepSizeStats,20);
    figure; subplot(8,1,1:5); plot(yInput,'b','LineWidth',2), hold on; plot(Step.StepFit,'r','LineWidth',2); 
    ylabel('y (nm)'); xlabel('Frame Number'); title(['On-axis step - ' FileInputName{ind}]);
    Stepind = 0;
    for i = 1:length(Step.StepSizeStats)
        Stepind = Stepind + Step.StepLengths(i);
        if Step.StepSizeStats(i) < 0
           text(Stepind+3,Step.StepFit(Stepind-1)+5,num2str(round(Step.StepSizeStats(i))),...
            'VerticalAlignment','middle',...
            'HorizontalAlignment','left',...
            'FontSize',8)
        elseif Step.StepSizeStats(i) > 0
            text(Stepind+3,Step.StepFit(Stepind-1)-5,num2str(round(Step.StepSizeStats(i))),...
            'VerticalAlignment','middle',...
            'HorizontalAlignment','left',...
            'FontSize',8)
        end
    end
    hold off;
    subplot(8,1,7:8); bar(xout,HistStep); ylabel('Count'); xlabel('Step-size (nm)'); title('Step-size histogram');
    cd(DataPath);
    save([FileInputName{ind} '-SICStepFinder.mat'],'Step');
    %% Analyze all stepping files
    for i = 1:length(FileIn)
        ind = i;
        fid=fopen(FileInput{ind});
        Input = textscan(fid,'%f%f','CommentStyle','##');
        yInput = Input{1};
        cd(CodePath);
        Step=SICstepFinder(yInput');
        figure; subplot(8,1,1:5); plot(yInput,'b','LineWidth',2), hold on; plot(Step.StepFit,'r','LineWidth',2); 
        ylabel('y (nm)'); xlabel('Frame Number'); title(['On-axis step - ' FileInputName{ind}]);
        Stepind = 0;
        for i = 1:length(Step.StepSizeStats)
            Stepind = Stepind + Step.StepLengths(i);
            if Step.StepSizeStats(i) < 0
               text(Stepind+3,Step.StepFit(Stepind-1)+5,num2str(round(Step.StepSizeStats(i))),...
                'VerticalAlignment','middle',...
                'HorizontalAlignment','left',...
                'FontSize',8)
            elseif Step.StepSizeStats(i) > 0
                text(Stepind+3,Step.StepFit(Stepind-1)-5,num2str(round(Step.StepSizeStats(i))),...
                'VerticalAlignment','middle',...
                'HorizontalAlignment','left',...
                'FontSize',8)
            end
        end
        [HistStep,xout] = hist(Step.StepSizeStats,20);
        subplot(8,1,7:8); bar(xout,HistStep); ylabel('Count'); xlabel('Step-size (nm)'); title('Step-size histogram');
        cd(DataPath);
        save([FileInputName{ind} '-SICStepFinder.mat'],'Step');
    end
    %% Analyze all on-axis stepping files and calculate FIONA Index
    for i = 1:length(FileIn)
        ind = i;
        fid=fopen(FileInput{ind});
        Input = textscan(fid,'%f%f','CommentStyle','##');
        yInput = Input{1};
        cd(CodePath);
        Step=SICstepFinder(yInput');
        [HistStep,xout] = hist(Step.StepSizeStats,100);
        cd(DataPath);
        figure; plot(yInput,'b','LineWidth',2), hold on; plot(Step.StepFit,'r','LineWidth',2); 
        Stepind = 0;
        for i = 1:length(Step.StepSizeStats)
            Stepind = Stepind + Step.StepLengths(i);
            if Step.StepSizeStats(i) < 0
               text(Stepind+3,Step.StepFit(Stepind-1)+5,num2str(round(Step.StepSizeStats(i))),...
                'VerticalAlignment','middle',...
                'HorizontalAlignment','left',...
                'FontSize',8)
            elseif Step.StepSizeStats(i) > 0
                text(Stepind+3,Step.StepFit(Stepind-1)-5,num2str(round(Step.StepSizeStats(i))),...
                'VerticalAlignment','middle',...
                'HorizontalAlignment','left',...
                'FontSize',8)
            end
        end
                
        % Calculate FIONA Index
        Stepsize = 8;          % Step size in nm
        FPS1 = length(yInput)*Stepsize/(max(yInput)-min(yInput));
        FPS2 = length(yInput)/length(Step.StepSizeStats);
        if FPS1 < FPS2; FPS = FPS1; FIndMethod = '(Distance)';
        else FPS = FPS2; FIndMethod = '(Steps)';
        end
        NPS = std(yInput-Step.StepFit')/Stepsize;
        FIONAIndex = FPS/(40*NPS^2+4.5*NPS+3.5);
        
        % Include these information in title
        title(['[FPS: ' sprintf('%1.1f',FPS) '] [NPS: ' sprintf('%1.1f',NPS) '] [FInd: ' sprintf('%1.1f',FIONAIndex) ']'],'fontweight','b','fontsize',12,'FontName','Palatino Linotype');
        
        % Save Image
        print('-djpeg','-r300',[FileInputName{ind} '-StepTrace.jpg']);
        
        hold off;
        figure; bar(xout,HistStep);
        print('-djpeg','-r300',[FileInputName{ind} '-Hist.jpg']);
        
        % Save *.mat files
        save(['FInd-' sprintf('%1.1f',FIONAIndex) '-' FIndMethod FileInputName{ind} '-SICStepFinder.mat'],'Step');
    end
    %% Analyze all off-axis stepping files (from Transformed file) and calculate FIONA Index
    for i = 1:length(FileIn)
        ind = i;
        fid=fopen(FileInput{ind});
        Input = textscan(fid,'%f%f%f','CommentStyle','##');
        yInput = Input{2};
        cd(CodePath);
        Step=SICstepFinder(yInput');
        [HistStep,xout] = hist(Step.StepSizeStats,100);
        cd(DataPath);
        figure; plot(yInput,'b','LineWidth',2), hold on; plot(Step.StepFit,'r','LineWidth',2); 
        Stepind = 0;
        for i = 1:length(Step.StepSizeStats)
            Stepind = Stepind + Step.StepLengths(i);
            if Step.StepSizeStats(i) < 0
               text(Stepind+3,Step.StepFit(Stepind-1)+5,num2str(round(Step.StepSizeStats(i))),...
                'VerticalAlignment','middle',...
                'HorizontalAlignment','left',...
                'FontSize',8)
            elseif Step.StepSizeStats(i) > 0
                text(Stepind+3,Step.StepFit(Stepind-1)-5,num2str(round(Step.StepSizeStats(i))),...
                'VerticalAlignment','middle',...
                'HorizontalAlignment','left',...
                'FontSize',8)
            end
        end
                
        % Calculate FIONA Index
        Stepsize = 8;          % Step size in nm
        FPS = length(yInput)*Stepsize/(max(yInput)-min(yInput));
        %FPS = length(yInput)/length(Step.StepSizeStats);
        NPS = std(yInput-Step.StepFit')/Stepsize;
        FIONAIndex = FPS/(40*NPS^2+4.5*NPS+3.5);
        
        % Include these information in title
        title(['[FPS: ' sprintf('%1.1f',FPS) '] [NPS: ' sprintf('%1.1f',NPS) '] [FInd: ' sprintf('%1.1f',FIONAIndex) ']'],'fontweight','b','fontsize',12,'FontName','Palatino Linotype');
        
        % Save Image
        print('-djpeg','-r300',[FileInputName{ind} '-StepTrace.jpg']);
        
        hold off;
        figure; bar(xout,HistStep);
        print('-djpeg','-r300',[FileInputName{ind} '-Hist.jpg']);
        
        % Save *.mat files
        save(['FInd-' sprintf('%1.1f',FIONAIndex) '-' FileInputName{ind} '-SICStepFinder.mat'],'Step');
    end
%% Draw step-size histogram after t-test analysis
    %% Search *.mat file in the folder (from SIC Stepfinder)
    clear
    % Input xThreshold to determine how wide plot range should be
    xThreshold = 0.99;
    
    % Search mat file and import data
    FileType = '.mat';
    FileIn = dir(['*' FileType]);
    FileInput = cell(length(FileIn),1);
    FileInputName = cell(length(FileIn),1);
    for i=1:length(FileIn)
        FileInput{i}=FileIn(i).name;
        FileInputName{i}=strrep(FileIn(i).name,FileType,'');
    end
    % Collect steps fom all *.mat files and save data
    CompiledDataPos = [];
    CompiledDataNeg = [];
    for i=1:length(FileIn)
        Step = load(FileInput{i}); Step = Step.Step;
        dataPos = Step.StepSizeStats(Step.StepSizeStats>0)';
        CompiledDataPos = [CompiledDataPos;dataPos];
        dataNeg = Step.StepSizeStats(Step.StepSizeStats<0)';
        CompiledDataNeg = [CompiledDataNeg;dataNeg];
    end
    % Add a line to prevent error when there is no data on negative side
    if isempty(CompiledDataNeg); CompiledDataNeg = -8; end
    %dlmwrite('CompiledPositiveSteps.txt',CompiledDataPos,'\t');
    %dlmwrite('CompiledNegativeSteps.txt',CompiledDataPos,'\t');
    
    % Get ranges
    binSpacing = 2;
    xHistPos = 0:binSpacing:ceil(max(CompiledDataPos));
    xHistNeg = floor(min(CompiledDataNeg)/binSpacing)*binSpacing:binSpacing:0;
    yHistPos=hist(CompiledDataPos,xHistPos);
    binSpace = xHistPos(2)-xHistPos(1);
    yHistNeg=hist(CompiledDataNeg,xHistNeg);
    CumSumNeg = cumsum(fliplr(yHistNeg)); 
    CumSumNeg = fliplr(CumSumNeg/CumSumNeg(end));
    NegInd = find(CumSumNeg>xThreshold,1,'last');
    CumSumPos = cumsum(yHistPos); CumSumPos = CumSumPos/CumSumPos(end);
    PosInd = find(CumSumPos>xThreshold,1);
    if abs(xHistNeg(NegInd))<xHistPos(PosInd)
        while abs(xHistNeg(NegInd))<xHistPos(PosInd) && NegInd > 1
            NegInd = NegInd - 1;
        end
    elseif abs(xHistNeg(NegInd))>xHistPos(PosInd)
        while abs(xHistNeg(NegInd))>xHistPos(PosInd) && PosInd <= length(xHistPos)
            PosInd = PosInd + 1;
        end
    end
    %% Gaussian Fitting (Two peaks at both -ve and +ve) 
    TextSize = 15;               % Size of text (for peak and stdev)
    xTickSpacing = 8;           % Specify the spacing for x-tick
    GaussColor = 'red';
    GraphTitle = 'K432 Step Size Histogram - Positive and Negative Step';
    %xHistPos = 0:binSpacing:ceil(max(CompiledDataPos));
    %xHistNeg = floor(min(CompiledDataNeg)/binSpacing)*binSpacing:binSpacing:0;
    yHistPos=hist(CompiledDataPos,xHistPos);
    yHistNeg=hist(CompiledDataNeg,xHistNeg);
    xHistCombined = [xHistNeg xHistPos(2:end)];
    yHistCombined=hist([CompiledDataPos; CompiledDataNeg], xHistCombined);
    %binSpace = xHistPos(2)-xHistPos(1);

    % Draw histogram
    h=bar(xHistCombined,yHistCombined,'hist');
    set(h,'facecolor',GaussColor,'LineWidth',1);
    sh=findall(gcf,'marker','*'); delete(sh); % delete the stars (*) that shows up in bar histogram
    xTickMin = floor(xHistNeg(NegInd)/xTickSpacing)*xTickSpacing;
    axGauss = gca; set(axGauss,'YColor','black','FontSize',TextSize,'XTick',xTickMin:xTickSpacing:xHistPos(PosInd));
    hold on;
    %bar(xHistNeg(NegInd:end),yHistNeg(NegInd:end),'r','LineWidth',1.5);

    center1 = 8;                    center2 = -8;
    std1 = 4;                       std2 = 4;
    peak1 = max(yHistPos);          peak2 = max(yHistNeg);
    
    center3 = 16;                   center4 = -16;
    std3 = 4;                       std4 = 4;
    peak3 = max(yHistPos)/4;        peak4 = max(yHistNeg)/4;
    %p = [peak1 center1 std1 peak2 center2 std2];
    %f = @(p,x)p(1)*exp(-((x-p(2))/p(3)).^2)+p(4)*exp(-((x-p(5))/p(6)).^2);
    %pfit = lsqcurvefit(f,p,xHist',yHist'); 

    %p1 = [peak1 center1 std1];
    %f = @(p,x)p(1)*exp(-((x-p(2))/(sqrt(2)*p(3))).^2);
    p1 = [peak1 center1 std1 peak3 center3 std3];
    f = @(p,x)p(1)*exp(-((x-p(2))/(sqrt(2)*p(3))).^2)+p(4)*exp(-((x-2*p(2))/(sqrt(2)*p(6))).^2);
    pfit1 = lsqcurvefit(f,p1,xHistPos',yHistPos'); 
    pfit1(5) = pfit1(2)*2;
    
    %p2 = [peak2 center2 std2];
    p2 = [peak2 center2 std2 peak4 center4 std4];
    pfit2 = lsqcurvefit(f,p2,xHistNeg',yHistNeg'); 
    pfit2(5) = pfit2(2)*2;

%     xfit = floor(min(xHistNeg)):0.5:ceil(max(xHistPos));
%     yfit1 = f(pfit1,xfit);
%     line(xfit,yfit1,'Color','black','LineWidth',3,'LineStyle','-', 'Parent',axGauss);
%     yfit2 = f(pfit2,xfit);
%     line(xfit,yfit2,'Color','black','LineWidth',3,'LineStyle','-', 'Parent',axGauss);
    f1gauss = @(p,x)p(1)*exp(-((x-p(2))/(sqrt(2)*p(3))).^2);
    xfit = floor(min(xHistNeg)):0.5:ceil(max(xHistPos));
    yfit1 = f1gauss(pfit1(1:3),xfit);
    yfit2 = f1gauss(pfit1(4:6),xfit);
    yfit3 = f1gauss(pfit2(1:3),xfit);
    yfit4 = f1gauss(pfit2(4:6),xfit);
    plot(xfit,yfit1,'--g','LineWidth',2);
    plot(xfit,yfit2,'--g','LineWidth',2);
    plot(xfit,yfit3,'--g','LineWidth',2);
    plot(xfit,yfit4,'--g','LineWidth',2);
    plot(xfit,yfit1+yfit2+yfit3+yfit4,'--b','LineWidth',3);
    
    ylabel('Count','fontweight','b','fontsize',TextSize,'FontName','Palatino Linotype');
    xlabel('Step size (nm)','fontweight','b','fontsize',12,'FontName','Palatino Linotype');
    title(GraphTitle,'fontweight','b','fontsize',TextSize,'FontName','Palatino Linotype');
    xlim([xHistNeg(NegInd)-center1/2 xHistPos(PosInd)+center1/2]);
    Text = str2mat(['Peak (-) = ' num2str(round(pfit2(2)*10)/10) ' ± ' num2str(round(pfit2(3))) ' nm'],['Peak (+) = ' num2str(round(pfit1(2)*10)/10) ' ± ' num2str(round(pfit1(3))) ' nm']);

    % Horizontal axis label (a textbox)
    PeakNeg = round(pfit2(2)*10)/10;
    PeakPos = round(pfit1(2)*10)/10; 
    nPos1 = integral(@(x)f1gauss(pfit1(1:3),x),0,max(abs(xHistPos)));
    semPos = pfit1(3)/sqrt(nPos1);
    nNeg1 = integral(@(x)f1gauss(pfit2(1:3),x),-max(abs(xHistNeg)),0);
    semNeg = pfit2(3)/sqrt(nNeg1);
    
    Negfract = sum(yHistNeg) * 100 / (sum(yHistNeg) + sum(yHistPos));
    Posfract = sum(yHistPos) * 100 / (sum(yHistNeg) + sum(yHistPos));
    
    TextNeg = str2mat([sprintf('%1.1f',PeakNeg) ' ± ' num2str(round((semNeg*10))/10) ' nm (' num2str(round(Negfract)) '%)']);
    text(PeakNeg-25,f(pfit2,PeakNeg)*1.2,TextNeg,'FontSize',TextSize,'fontweight','b')
    TextPos = str2mat([sprintf('%1.1f',PeakPos)  ' ± ' num2str(round((semPos*10))/10) ' nm (' num2str(round(Posfract)) '%)']);
    text(PeakPos-3,f(pfit1,PeakPos)*1.05,TextPos,'FontSize',TextSize,'fontweight','b')

    print('-djpeg','-r300',[GraphTitle '.png']);




%f = @(x,xdata)x(1)*xdata.^2+x(2)*sin(xdata);
%x = lsqcurvefit(f,x0,xdata,ydata);
%a1*exp(-((x-b1)/c1)^2)+a2*exp(-((x-b2)/c2)^2;
    %% Gaussian Fitting (Two peaks at +ve)
    TextSize = 15;               % Size of text (for peak and stdev)
    xTickSpacing = 8;           % Specify the spacing for x-tick
    binSpacing = 2;
    GraphTitle = 'K432 Step Size Histogram - Two Gaussian at Positive Steps';
    xHistPos = 0:binSpacing:ceil(max(CompiledDataPos));
    xHistNeg = floor(min(CompiledDataNeg)/binSpacing)*binSpacing:binSpacing:0;
    yHistPos=hist(CompiledDataPos,xHistPos);

    xHistPos = xHistPos(yHistPos~=0); yHistPos = yHistPos(yHistPos~=0);

    figure;
    bar(xHistPos,yHistPos,'r','LineWidth',1.5);
    axGauss = gca; set(axGauss,'YColor','black','FontSize',TextSize,'XTick',0:xTickSpacing:xHistPos(PosInd));
    hold on;

    center1 = 8;                    center2 = 16;
    std1 = 4;                      std2 = 4;
    peak1 = max(yHistPos);          peak2 = max(yHistPos)/4;
    %p = [peak1 center1 std1 peak2 center2 std2];
    p = [peak1 center1 std1 peak2 2*center1 std2];    % p for second step size twice as big as first step size
    %f = @(p,x)p(1)*exp(-((x-p(2))/(sqrt(2)*p(3))).^2)+p(4)*exp(-((x-p(5))/(sqrt(2)*p(6))).^2);
    f = @(p,x)p(1)*exp(-((x-p(2))/(sqrt(2)*p(3))).^2)+p(4)*exp(-((x-(2*p(2)))/(sqrt(2)*p(6))).^2);  % f for second step size twice as big as first step size
    %f=@(p,x)(p(1)/(p(3)*sqrt(pi/2)))*sqrt(x/p(2)).*exp(-2*((x-p(2))/p(3)).^2)+...
    %    (p(4)/(p(6)*sqrt(pi/2)))*sqrt(x/p(5)).*exp(-2*((x-p(5))/p(6)).^2);
    %f=@(p,x)(p(1)/(p(3)*sqrt(pi/2)))*sqrt(x/p(2)).*exp(-2*((x-p(2))/p(3)).^2)+...
    %    (p(4)/(p(6)*sqrt(pi/2)))*sqrt(x/(2*p(2))).*exp(-2*((x-(2*p(2)))/p(6)).^2);
    pfit1 = lsqcurvefit(f,p,xHistPos',yHistPos'); 
    pfit1(5) = pfit1(2)*2;

    %p1 = [peak1 center1 std1];
    %f = @(p,x)p(1)*exp(-((x-p(2))/(sqrt(2)*p(3))).^2);
    %pfit1 = lsqcurvefit(f,p1,xHistPos',yHistPos'); 

    %p2 = [peak2 center2 std2];
    %pfit2 = lsqcurvefit(f,p2,xHistNeg',yHistNeg'); 
    
    f1gauss = @(p,x)p(1)*exp(-((x-p(2))/(sqrt(2)*p(3))).^2);
    xfit = floor(min(xHistNeg)):0.5:ceil(max(xHistPos));
    yfit1 = f1gauss(pfit1(1:3),xfit);
    yfit2 = f1gauss(pfit1(4:6),xfit);
    plot(xfit,yfit1,'--g','LineWidth',2);
    plot(xfit,yfit2,'--g','LineWidth',2);
    plot(xfit,yfit1+yfit2,'--b','LineWidth',3);
    %yfitAll = f(pfit1,xfit);
    %plot(xfit,yfitAll,'black','LineWidth',2);
    %yfit2 = f(pfit2,xfit);
    %plot(xfit,yfit2,'blue','LineWidth',2);
    
    % Horizontal axis label (a textbox)
    PeakNeg = round(pfit1(2)*10)/10;
    %nNeg = integral(@(x)f1gauss(pfit1(1:3),x),0,max(abs(xHistNeg)));
    nNeg = integral(@(x)f1gauss(pfit1(1:3),x),0,max(abs(xHistPos)));
    semNeg = pfit1(3)/sqrt(nNeg);
    PeakPos = round(pfit1(5)*10)/10; 
    nPos = integral(@(x)f1gauss(pfit1(4:6),x),0,max(xHistPos));
    semPos = pfit1(6)/sqrt(nPos);
    
    Negfract = nNeg * 100 / (nNeg + nPos);
    Posfract = nPos * 100 / (nNeg + nPos);

    TextNeg = str2mat([sprintf('%1.1f',PeakNeg) ' ± ' num2str(round((semNeg*10))/10) ' nm (' num2str(round(Negfract)) '%)']);
    text(PeakNeg+pfit1(3),f(pfit1,PeakNeg)*1,TextNeg,'FontSize',TextSize,'fontweight','b')
    TextPos = str2mat([sprintf('%1.1f',PeakPos)  ' ± ' num2str(round((semPos*10))/10) ' nm (' num2str(round(Posfract)) '%)']);
    text(PeakPos+pfit1(6),f(pfit1,PeakPos)*1,TextPos,'FontSize',TextSize,'fontweight','b')

    ylabel('Frequency','fontweight','b','fontsize',12,'FontName','Palatino Linotype');
    xlabel('Step size (nm)','fontweight','b','fontsize',12,'FontName','Palatino Linotype');
    title([GraphTitle],'fontweight','b','fontsize',12,'FontName','Palatino Linotype');
    xlim([0 xHistPos(PosInd)]);

    %print('-djpeg','-r300',[GraphTitle '.jpg']);
    
    % Plot residual
    yfit1 = f1gauss(pfit1(1:3),xHistPos);
    yfit2 = f1gauss(pfit1(4:6),xHistPos);
    bar(xHistPos,yfit1+yfit2-yHistPos);
    print('-djpeg','-r300',[GraphTitle '.jpg']);    
%% Draw dwell-time histogram after t-test analysis   
    %% Search *.mat file in the folder (from SIC Stepfinder in T-test folder)
    clear
    FileType = '.mat';
    FileIn = dir(['*' FileType]);
    FileInput = cell(length(FileIn),1);
    FileInputName = cell(length(FileIn),1);
    for i=1:length(FileIn)
        FileInput{i}=FileIn(i).name;
        FileInputName{i}=strrep(FileIn(i).name,FileType,'');
        FilePrefix(i)=strrep(FileInputName(i),'-SICStepFinder','');
    end
    
    % Collect dwell time data from all files
    CompiledDwell = [];
    for i=1:length(FileIn)
        fid=fopen([FilePrefix{i} '.txt']);
        Input = textscan(fid,'%f%f','CommentStyle','##');
        tInput = Input{2};
        Step = load(FileInput{i}); Step = Step.Step;
        StepFit = Step.StepFit'; StepFit = StepFit(2:end) - StepFit(1:end-1);
        dwell = tInput(find(StepFit~=0)); dwell = dwell(2:end) - dwell(1:end-1);
        CompiledDwell = [CompiledDwell;dwell];
    end
    %% Plot and fit dwell histogram (hidden step) (Combine -ve and +ve)
    ExposureTime = 0.05;            % Exposure time in second
    binSpacing = 4;
    GraphTitle = 'Kinesin Dwell Histogram';
    xHistDwell = 0:binSpacing:ceil(max(CompiledDwell));
    yHistDwell=hist(CompiledDwell,xHistDwell);

    binSpace = xHistDwell(2)-xHistDwell(1);
    xHistDwell = xHistDwell(yHistDwell~=0)*ExposureTime; yHistDwell = yHistDwell(yHistDwell~=0);

    figure;
    bar(xHistDwell,yHistDwell,'r','LineWidth',1.5);
    hold on;
    
    k = 0.05;
    A = 30;
    p = [A k];
    f = @(p,x)p(1)*x.*exp(-p(2)*x);
    pfit = lsqcurvefit(f,p,xHistDwell',yHistDwell'); 

    xfit = 0:0.05:ceil(max(xHistDwell));
    yfit = f(pfit,xfit);
    plot(xfit,yfit,'blue','LineWidth',2);


    ylabel('Frequency','fontweight','b','fontsize',20,'FontName','Palatino Linotype');
    xlabel('Dwell time (second)','fontweight','b','fontsize',20,'FontName','Palatino Linotype');
    title(GraphTitle,'fontweight','b','fontsize',20,'FontName','Palatino Linotype');
    set(gca,'fontsize',20);
    xlim([0 5]);
    Text = str2mat(['k = ' num2str(pfit(2)) ' s^-^1']);

    % Horizontal axis label (a textbox)
    annotation('textbox',[0.5 0.8 0.8 0.05],'String',Text,...
    'FontWeight','bold','FontSize',20,'FontName','Palatino Linotype',...
    'FontAngle','italic','FitBoxToText','off','LineStyle','none');

    print('-djpeg','-r300',[GraphTitle '.jpg']);
    %% Plot and fit dwell histogram (exponential) (Combine -ve and +ve)
    ExposureTime = 0.05;            % Exposure time in second
    binSpacing = 4;
    GraphTitle = 'Kinesin Dwell Histogram';
    xHistDwell = 0:binSpacing:ceil(max(CompiledDwell));
    yHistDwell=hist(CompiledDwell,xHistDwell);

    binSpace = xHistDwell(2)-xHistDwell(1);
    xHistDwell = xHistDwell(yHistDwell~=0)*ExposureTime; yHistDwell = yHistDwell(yHistDwell~=0);

    figure;
    bar(xHistDwell,yHistDwell,'r','LineWidth',1.5);
    hold on;

    k = 0.05;
    A = 30;
    p = [A k];
    f = @(p,x)p(1)*exp(-p(2)*x);
    pfit = lsqcurvefit(f,p,xHistDwell',yHistDwell'); 

    xfit = 0:0.05:ceil(max(xHistDwell));
    yfit = f(pfit,xfit);
    plot(xfit,yfit,'blue','LineWidth',2);


    ylabel('Frequency','fontweight','b','fontsize',20,'FontName','Palatino Linotype');
    xlabel('Dwell time (second)','fontweight','b','fontsize',20,'FontName','Palatino Linotype');
    title(GraphTitle,'fontweight','b','fontsize',20,'FontName','Palatino Linotype');
    set(gca,'fontsize',20);
    xlim([0 5]);
    Text = str2mat(['k = ' num2str(pfit(2)) ' s^-^1']);

    % Horizontal axis label (a textbox)
    annotation('textbox',[0.5 0.8 0.8 0.05],'String',Text,...
    'FontWeight','bold','FontSize',20,'FontName','Palatino Linotype',...
    'FontAngle','italic','FitBoxToText','off','LineStyle','none');

    print('-djpeg','-r300',[GraphTitle '.jpg']);
%% Calculate velocity and run length from T-test folder
ExposureTime = 0.1;     % Exposure time in unit of second
FileType = '.txt';
FileIn = dir(['*' FileType]);
FileInput = cell(length(FileIn),1);
FileInputName = cell(length(FileIn),1);
for ind=1:length(FileIn)
    FileInput{ind}=FileIn(ind).name;
    FileInputName{ind}=strrep(FileIn(ind).name,FileType,'');
end
% Initiate velocity
Vel = zeros(length(FileIn),1);
RunLength = zeros(length(FileIn),1);

for ind = 1:length(FileIn)
    fid=fopen(FileInput{ind});
    Input = textscan(fid,'%f%f','CommentStyle','##');
    yInput = Input{1};
    figure; plot(yInput); hold on;
    plot([0 length(yInput)], [min(yInput) max(yInput)],'r','LineWidth',3);
    Vel(ind) = (max(yInput)-min(yInput))/(length(yInput)*ExposureTime);     % Velocity in nm/s
    RunLength(ind) = (max(yInput)-min(yInput));     % Run Length in nm
    title(['[Velocity = ' num2str(Vel(ind)) ' nm/s] [Run Length = ' num2str(RunLength(ind)) ' nm]']);
end

figure; subplot(2,1,1); hist(Vel,5); MeanVel = mean(Vel);
title(['[Mean Velocity = ' num2str(round(MeanVel*100)/100) ' nm/s] [Standard Deviation = ' num2str(round(std(Vel)*100)/100) ' nm/s]']);
subplot(2,1,2); hist(RunLength,5); MeanRunLength = mean(RunLength);
title(['[Mean Run Length = ' num2str(round(MeanRunLength*10)/10) ' nm] [Standard Deviation = ' num2str(round(std(RunLength)*10)/10) ' nm]']);





%% More Codes
    %% Two modified gaussian peaks at positive steps
    NegShift = 7;
    PosShift = 5;
    TextSize = 15;               % Size of text (for peak and stdev)
    xTickSpacing = 8;           % Specify the spacing for x-tick
    binSpacing = 2;
    GraphTitle = 'DDB Step Size Histogram - Two Modified Gaussian at Positive Steps';
    xHistPos = 0:binSpacing:ceil(max(CompiledDataPos));
    xHistNeg = floor(min(CompiledDataNeg)/binSpacing)*binSpacing:binSpacing:0;
    yHistPos=hist(CompiledDataPos,xHistPos);

    binSpace = xHistPos(2)-xHistPos(1);
    yHistNeg=hist(CompiledDataNeg,xHistNeg);

    %xHistNeg = xHist(xHist<0); yHistNeg = yHist(xHist<0);
    xHistNeg = xHistNeg(yHistNeg~=0); yHistNeg = yHistNeg(yHistNeg~=0);
    %xHistPos = xHist(xHist>=0); yHistPos = yHist(xHist>=0);
    xHistPos = xHistPos(yHistPos~=0); yHistPos = yHistPos(yHistPos~=0);

    figure;
    bar(xHistPos,yHistPos,'r','LineWidth',1.5);
    hold on;
    %bar(xHistNeg,yHistNeg,'r','LineWidth',1.5);

    center1 = 8;            center2 = 16;
    std1 = 10;               std2 = 10;
    peak1 = 170;            peak2 = 20;
    p = [peak1 center1 std1 peak2 center2 std2];
    %f=@(p,x)(p(1)/(p(3)*sqrt(pi/2)))*sqrt(x/p(2)).*exp(-2*((x-p(2))/p(3)).^2)+...
    %    (p(4)/(p(6)*sqrt(pi/2)))*sqrt(x/p(5)).*exp(-2*((x-p(5))/p(6)).^2);
    f=@(p,x)(p(1)/(p(3)*sqrt(pi/2)))*sqrt(x/p(2)).*exp(-2*((x-p(2))/p(3)).^2)+...
        (p(4)/(p(6)*sqrt(pi/2)))*sqrt(x/(2*p(2))).*exp(-2*((x-(2*p(2)))/p(6)).^2);
    pfit1 = lsqcurvefit(f,p,xHistPos',yHistPos'); 
    pfit1(5) = pfit1(2)*2;
    %p1 = [peak1 center1 std1];
    %f = @(p,x)p(1)*exp(-((x-p(2))/(sqrt(2)*p(3))).^2);
    %pfit1 = lsqcurvefit(f,p1,xHistPos',yHistPos'); 

    %p2 = [peak2 center2 std2];
    %pfit2 = lsqcurvefit(f,p2,xHistNeg',yHistNeg'); 

    f1gauss = @(p,x)(p(1)/(p(3)*sqrt(pi/2)))*sqrt(x/p(2)).*exp(-2*((x-p(2))/p(3)).^2);
    xfit = floor(min(xHistNeg)):0.5:ceil(max(xHistPos));
    yfit1 = f1gauss(pfit1(1:3),xfit);
    yfit2 = f1gauss(pfit1(4:6),xfit);
    plot(xfit,yfit1,'--g','LineWidth',2);
    plot(xfit,yfit2,'--g','LineWidth',2);
    plot(xfit,yfit1+yfit2,'--b','LineWidth',3);
    %yfitAll = f(pfit1,xfit);
    %plot(xfit,yfitAll,'black','LineWidth',2);
    %yfit2 = f(pfit2,xfit);
    %plot(xfit,yfit2,'blue','LineWidth',2);

    % Horizontal axis label (a textbox)
    PeakNeg = round(pfit1(2)*10)/10;
    nNeg = integral(@(x)f1gauss(pfit1(1:3),x),0,max(abs(xHistNeg)));
    semNeg = pfit1(3)/sqrt(nNeg);
    TextNeg = str2mat([sprintf('%1.1f',PeakNeg) ' ± ' num2str(round((semNeg*10))/10) ' nm']);
    text(PeakNeg+NegShift,f(pfit1,PeakNeg)*1,TextNeg,'FontSize',TextSize,'fontweight','b')
    PeakPos = round(pfit1(5)*10)/10; 
    nPos = integral(@(x)f1gauss(pfit1,x),0,max(xHistPos));
    semPos = pfit1(6)/sqrt(nPos);
    TextPos = str2mat([sprintf('%1.1f',PeakPos)  ' ± ' num2str(round((semPos*10))/10) ' nm']);
    text(PeakPos+PosShift,f(pfit1,PeakPos)*1,TextPos,'FontSize',TextSize,'fontweight','b')

    ylabel('Frequency','fontweight','b','fontsize',12,'FontName','Palatino Linotype');
    xlabel('Step size (nm)','fontweight','b','fontsize',12,'FontName','Palatino Linotype');
    title([GraphTitle],'fontweight','b','fontsize',12,'FontName','Palatino Linotype');
    xlim([0 xHistPos(PosInd)]);

    % Plot residual
    yfit1 = f1gauss(pfit1(1:3),xHistPos);
    yfit2 = f1gauss(pfit1(4:6),xHistPos);
    bar(xHistPos,yfit1+yfit2-yHistPos);
    
    % Save image
    print('-djpeg','-r300',[GraphTitle '.jpg']);
    %% Maximum likelihood with automatic search for Peak No (Combine -ve and +ve)
    binSpacing = 1;
    xThreshold = 0.985;
    AnnotateXLocation = 0.45;    % Range from 0 to 1. Determine horizontal location of annotation
    AnnotationSize = 15;
    GraphTitle = 'Kinesin Step Size Histogram';
    CompiledData = [CompiledDataPos;-(CompiledDataNeg)];
    xHist = 0:binSpacing:ceil(max(CompiledData));
    yHist=hist(CompiledData,xHist);
    
    % Find out variables to make the graph look full
    CumSumPos = cumsum(yHist); CumSumPos = CumSumPos/CumSumPos(end);
    PosInd = find(CumSumPos>xThreshold,1);

    figure;
    bar(xHist(1:PosInd),yHist(1:PosInd),'r','LineWidth',1.5);
    hold on;

    options = statset('Display','final');
    Length = length(CompiledData);
    xfit = xHist'; PeakNo = 1; RSquarePre = 0;
    obj = gmdistribution.fit(CompiledData,1,'Options',options);
    yfit = (pdf(obj,xfit)*Length)';
    resnorm = sum((yHist-yfit).^2);
    SumSquaresTotal = sum((yHist-mean(yHist)).^2);
    RSquare = 1 - (resnorm/SumSquaresTotal);
    BICPre = obj.BIC;

    PeakNo = 2;
    obj = gmdistribution.fit(Step.StepSizeStats',PeakNo,'Options',options);
    BIC = obj.BIC;
    while BICPre > BIC
        PeakNo = PeakNo + 1;
        obj = gmdistribution.fit(Step.StepSizeStats',PeakNo,'Options',options);
        BICPre = BIC;
        BIC = obj.BIC;
    end

    obj = gmdistribution.fit(CompiledData,PeakNo-1,'Options',options);
    yfit = (pdf(obj,xfit)*Length)';

    plot(xfit,yfit*binSpacing,'--b','LineWidth',3);
    ylabel('Frequency','fontweight','b','fontsize',12,'FontName','Palatino Linotype');
    xlabel('Step size (nm)','fontweight','b','fontsize',12,'FontName','Palatino Linotype');
    title([GraphTitle ' - Max Likelihood Fit (all)'],'fontweight','b','fontsize',12,'FontName','Palatino Linotype');
    [Peaks, IX] = sort(round(obj.mu*10)/10);
    f1gauss = @(p,x)p(1)*exp(-((x-p(2))/(sqrt(2)*p(3))).^2);
    Std = sqrt(obj.Sigma(1));
    NoOfPeaks = PeakNo-1;
    if NoOfPeaks > 1;
        for iter = 2:NoOfPeaks; Std = [Std sqrt(obj.Sigma(iter))]; end
    end
    Amplitude = obj.PComponents*Length./(Std*sqrt(2*pi))*binSpacing;
    nGauss = Amplitude.*Std*sqrt(2*pi);
    sem = round(Std./sqrt(nGauss)*10)/10;
    yfitgauss = zeros(length(xfit),NoOfPeaks);
    
    Text = str2mat(['Peak 1 = ' num2str(Peaks(1)) ' ± ' num2str(sem(IX(1))) ' nm']);
    yfitgauss(:,IX(1)) = f1gauss([Amplitude(IX(1)) obj.mu(IX(1)) Std(IX(1))],xfit);
    plot(xfit,yfitgauss(:,IX(1)),'--g','LineWidth',2);
    if NoOfPeaks > 1
        for iter = 2:NoOfPeaks
            Text = str2mat(Text,['Peak ' num2str(iter) ' = ' num2str(Peaks(iter)) ' ± ' num2str(sem(IX(iter))) ' nm']);
            yfitgauss(:,IX(iter)) = f1gauss([Amplitude(IX(iter)) obj.mu(IX(iter)) Std(IX(iter))],xfit);
            plot(xfit,yfitgauss(:,IX(iter)),'--g','LineWidth',2);
        end
    end
    
    % Horizontal axis label (a textbox)
    annotation('textbox',[AnnotateXLocation 0.8 0.8 0.05],'String',Text,...
    'FontWeight','bold','FontSize',AnnotationSize,'FontName','Palatino Linotype',...
    'FontAngle','italic','FitBoxToText','off','LineStyle','none');
    xlim([0 xHist(PosInd)]);

    print('-djpeg','-r300',[GraphTitle '-MLE-All.jpg']);    
    %% Maximum likelihood with automatic search for Peak No (+ve only)
    binSpacing = 1;
    xThreshold = 0.985;
    AnnotateXLocation = 0.45;    % Range from 0 to 1. Determine horizontal location of annotation
    AnnotationSize = 15;
    GraphTitle = 'Kinesin Step Size Histogram';
    CompiledData = CompiledDataPos;
    xHist = 0:binSpacing:ceil(max(CompiledData));
    yHist=hist(CompiledData,xHist);
    
    % Find out variables to make the graph look full
    CumSumPos = cumsum(yHist); CumSumPos = CumSumPos/CumSumPos(end);
    PosInd = find(CumSumPos>xThreshold,1);

    figure;
    bar(xHist(1:PosInd),yHist(1:PosInd),'r','LineWidth',1.5);
    hold on;

    options = statset('Display','final');
    Length = length(CompiledData);
    xfit = xHist'; PeakNo = 1; RSquarePre = 0;
    obj = gmdistribution.fit(CompiledData,1,'Options',options);
    yfit = (pdf(obj,xfit)*Length)';
    resnorm = sum((yHist-yfit).^2);
    SumSquaresTotal = sum((yHist-mean(yHist)).^2);
    RSquare = 1 - (resnorm/SumSquaresTotal);
    BICPre = obj.BIC;

    PeakNo = 2;
    obj = gmdistribution.fit(Step.StepSizeStats',PeakNo,'Options',options);
    BIC = obj.BIC;
    while BICPre > BIC
        PeakNo = PeakNo + 1;
        obj = gmdistribution.fit(Step.StepSizeStats',PeakNo,'Options',options);
        BICPre = BIC;
        BIC = obj.BIC;
    end

    obj = gmdistribution.fit(CompiledData,PeakNo-1,'Options',options);
    yfit = (pdf(obj,xfit)*Length)';

    plot(xfit,yfit*binSpacing,'--b','LineWidth',3);
    ylabel('Frequency','fontweight','b','fontsize',12,'FontName','Palatino Linotype');
    xlabel('Step size (nm)','fontweight','b','fontsize',12,'FontName','Palatino Linotype');
    title([GraphTitle ' - Max Likelihood Fit (all)'],'fontweight','b','fontsize',12,'FontName','Palatino Linotype');
    [Peaks, IX] = sort(round(obj.mu*10)/10);
    f1gauss = @(p,x)p(1)*exp(-((x-p(2))/(sqrt(2)*p(3))).^2);
    Std = sqrt(obj.Sigma(1));
    NoOfPeaks = PeakNo-1;
    if NoOfPeaks > 1;
        for iter = 2:NoOfPeaks; Std = [Std sqrt(obj.Sigma(iter))]; end
    end
    Amplitude = obj.PComponents*Length./(Std*sqrt(2*pi))*binSpacing;
    nGauss = Amplitude.*Std*sqrt(2*pi);
    sem = round(Std./sqrt(nGauss)*10)/10;
    yfitgauss = zeros(length(xfit),NoOfPeaks);
    
    Text = str2mat(['Peak 1 = ' num2str(Peaks(1)) ' ± ' num2str(sem(IX(1))) ' nm']);
    yfitgauss(:,IX(1)) = f1gauss([Amplitude(IX(1)) obj.mu(IX(1)) Std(IX(1))],xfit);
    plot(xfit,yfitgauss(:,IX(1)),'--g','LineWidth',2);
    if NoOfPeaks > 1
        for iter = 2:NoOfPeaks
            Text = str2mat(Text,['Peak ' num2str(iter) ' = ' num2str(Peaks(iter)) ' ± ' num2str(sem(IX(iter))) ' nm']);
            yfitgauss(:,IX(iter)) = f1gauss([Amplitude(IX(iter)) obj.mu(IX(iter)) Std(IX(iter))],xfit);
            plot(xfit,yfitgauss(:,IX(iter)),'--g','LineWidth',2);
        end
    end
    
    % Horizontal axis label (a textbox)
    annotation('textbox',[AnnotateXLocation 0.8 0.8 0.05],'String',Text,...
    'FontWeight','bold','FontSize',AnnotationSize,'FontName','Palatino Linotype',...
    'FontAngle','italic','FitBoxToText','off','LineStyle','none');
    xlim([0 xHist(PosInd)]);

    print('-djpeg','-r300',[GraphTitle '-MLE-All.jpg']);    
    %% Maximum likelihood with manual defintion for Peak No (Combine -ve and +ve)
    binSpacing = 1;
    NoOfPeaks = 3;
    xThreshold = 0.985;
    AnnotateXLocation = 0.45;    % Range from 0 to 1. Determine horizontal location of annotation
    AnnotationSize = 15;
    GraphTitle = 'Kinesin Step Size Histogram';
    CompiledData = [CompiledDataPos;-(CompiledDataNeg)];
    xHist = 0:binSpacing:ceil(max(CompiledData));
    yHist=hist(CompiledData,xHist);
    
    % Find out variables to make the graph look full
    CumSumPos = cumsum(yHist); CumSumPos = CumSumPos/CumSumPos(end);
    PosInd = find(CumSumPos>xThreshold,1);

    figure;
    bar(xHist(1:PosInd),yHist(1:PosInd),'r','LineWidth',1.5);
    hold on;

    options = statset('Display','final');
    Length = length(CompiledData);
    xfit = xHist'; 

    obj = gmdistribution.fit(CompiledData,NoOfPeaks,'Options',options);
    yfit = (pdf(obj,xfit)*Length)';
    %yfit2 = (pdf(obj,(0:100)'))';
    
    plot(xfit,yfit*binSpacing,'--b','LineWidth',3);
    ylabel('Frequency','fontweight','b','fontsize',12,'FontName','Palatino Linotype');
    xlabel('Step size (nm)','fontweight','b','fontsize',12,'FontName','Palatino Linotype');
    title([GraphTitle ' - Max Likelihood Fit (all)'],'fontweight','b','fontsize',12,'FontName','Palatino Linotype');
    [Peaks, IX] = sort(round(obj.mu*10)/10);
    %Stdev = round(obj.Sigma);
    f1gauss = @(p,x)p(1)*exp(-((x-p(2))/(sqrt(2)*p(3))).^2);
    Std = sqrt(obj.Sigma(1));
    if NoOfPeaks > 1;
        for iter = 2:NoOfPeaks; Std = [Std sqrt(obj.Sigma(iter))]; end
    end
    Amplitude = obj.PComponents*Length./(Std*sqrt(2*pi))*binSpacing;
    nGauss = Amplitude.*Std*sqrt(2*pi);
    sem = round(Std./sqrt(nGauss)*10)/10;
    yfitgauss = zeros(length(xfit),NoOfPeaks);
    
    Text = str2mat(['Peak 1 = ' num2str(Peaks(1)) ' ± ' num2str(sem(IX(1))) ' nm']);
    yfitgauss(:,IX(1)) = f1gauss([Amplitude(IX(1)) obj.mu(IX(1)) Std(IX(1))],xfit);
    plot(xfit,yfitgauss(:,IX(1)),'--g','LineWidth',2);
    if NoOfPeaks > 1
        for iter = 2:NoOfPeaks
            Text = str2mat(Text,['Peak ' num2str(iter) ' = ' num2str(Peaks(iter)) ' ± ' num2str(sem(IX(iter))) ' nm']);
            yfitgauss(:,IX(iter)) = f1gauss([Amplitude(IX(iter)) obj.mu(IX(iter)) Std(IX(iter))],xfit);
            plot(xfit,yfitgauss(:,IX(iter)),'--g','LineWidth',2);
        end
    end
    
    % Horizontal axis label (a textbox)
    annotation('textbox',[AnnotateXLocation 0.8 0.8 0.05],'String',Text,...
    'FontWeight','bold','FontSize',AnnotationSize,'FontName','Palatino Linotype',...
    'FontAngle','italic','FitBoxToText','off','LineStyle','none');
    xlim([0 xHist(PosInd)]);

    print('-djpeg','-r300',[GraphTitle '-MLE-All.jpg']); 
    %% Maximum likelihood with manual defintion for Peak No (+ve only)
    binSpacing = 2;
    NoOfPeaks = 3;
    xThreshold = 0.985;
    AnnotateXLocation = 0.45;    % Range from 0 to 1. Determine horizontal location of annotation
    AnnotationSize = 15;
    GraphTitle = 'Kinesin Step Size Histogram';
    CompiledData = CompiledDataPos;
    xHist = 0:binSpacing:ceil(max(CompiledData));
    yHist=hist(CompiledData,xHist);
    
    % Find out variables to make the graph look full
    CumSumPos = cumsum(yHist); CumSumPos = CumSumPos/CumSumPos(end);
    PosInd = find(CumSumPos>xThreshold,1);

    figure;
    bar(xHist(1:PosInd),yHist(1:PosInd),'r','LineWidth',1.5);
    hold on;

    options = statset('Display','final');
    Length = length(CompiledData);
    xfit = xHist'; 

    obj = gmdistribution.fit(CompiledData,NoOfPeaks,'Options',options);
    yfit = (pdf(obj,xfit)*Length)';
    %yfit2 = (pdf(obj,(0:100)'))';
    
    plot(xfit,yfit*binSpacing,'--b','LineWidth',3);
    ylabel('Frequency','fontweight','b','fontsize',12,'FontName','Palatino Linotype');
    xlabel('Step size (nm)','fontweight','b','fontsize',12,'FontName','Palatino Linotype');
    title([GraphTitle ' - Max Likelihood Fit (all)'],'fontweight','b','fontsize',12,'FontName','Palatino Linotype');
    [Peaks, IX] = sort(round(obj.mu*10)/10);
    %Stdev = round(obj.Sigma);
    f1gauss = @(p,x)p(1)*exp(-((x-p(2))/(sqrt(2)*p(3))).^2);
    Std = sqrt(obj.Sigma(1));
    if NoOfPeaks > 1;
        for iter = 2:NoOfPeaks; Std = [Std sqrt(obj.Sigma(iter))]; end
    end
    Amplitude = obj.PComponents*Length./(Std*sqrt(2*pi))*binSpacing;
    nGauss = Amplitude.*Std*sqrt(2*pi);
    sem = round(Std./sqrt(nGauss)*10)/10;
    yfitgauss = zeros(length(xfit),NoOfPeaks);
    
    Text = str2mat(['Peak 1 = ' num2str(Peaks(1)) ' ± ' num2str(sem(IX(1))) ' nm']);
    yfitgauss(:,IX(1)) = f1gauss([Amplitude(IX(1)) obj.mu(IX(1)) Std(IX(1))],xfit);
    plot(xfit,yfitgauss(:,IX(1)),'--g','LineWidth',2);
    if NoOfPeaks > 1
        for iter = 2:NoOfPeaks
            Text = str2mat(Text,['Peak ' num2str(iter) ' = ' num2str(Peaks(iter)) ' ± ' num2str(sem(IX(iter))) ' nm']);
            yfitgauss(:,IX(iter)) = f1gauss([Amplitude(IX(iter)) obj.mu(IX(iter)) Std(IX(iter))],xfit);
            plot(xfit,yfitgauss(:,IX(iter)),'--g','LineWidth',2);
        end
    end
    
    % Horizontal axis label (a textbox)
    annotation('textbox',[AnnotateXLocation 0.8 0.8 0.05],'String',Text,...
    'FontWeight','bold','FontSize',AnnotationSize,'FontName','Palatino Linotype',...
    'FontAngle','italic','FitBoxToText','off','LineStyle','none');
    xlim([0 xHist(PosInd)]);

    print('-djpeg','-r300',[GraphTitle '-MLE-All.jpg']); 
    %% Gaussian Fitting (Separate -ve and +ve) 
    TextSize = 20;               % Size of text (for peak and stdev)
    xTickSpacing = 8;           % Specify the spacing for x-tick
    GaussColor = 'blue';
    GraphTitle = 'DDB Step Size Histogram';
    %xHistPos = 0:binSpacing:ceil(max(CompiledDataPos));
    %xHistNeg = floor(min(CompiledDataNeg)/binSpacing)*binSpacing:binSpacing:0;
    yHistPos=hist(CompiledDataPos,xHistPos);
    yHistNeg=hist(CompiledDataNeg,xHistNeg);
    xHistCombined = [xHistNeg xHistPos(2:end)];
    yHistCombined=hist([CompiledDataPos; CompiledDataNeg], xHistCombined);
    %binSpace = xHistPos(2)-xHistPos(1);

    % Draw histogram
    h=bar(xHistCombined,yHistCombined,'hist');
    set(h,'facecolor',GaussColor,'LineWidth',1.5);
    sh=findall(gcf,'marker','*'); delete(sh); % delete the stars (*) that shows up in bar histogram
    xTickMin = floor(xHistNeg(NegInd)/xTickSpacing)*xTickSpacing;
    axGauss = gca; set(axGauss,'YColor','black','FontSize',TextSize,'XTick',xTickMin:xTickSpacing:xHistPos(PosInd));
    hold on;
    %bar(xHistNeg(NegInd:end),yHistNeg(NegInd:end),'r','LineWidth',1.5);

    center1 = 8;            center2 = -8;
    std1 = 10;               std2 = 10;
    peak1 = 100;            peak2 = 10;
    %p = [peak1 center1 std1 peak2 center2 std2];
    %f = @(p,x)p(1)*exp(-((x-p(2))/p(3)).^2)+p(4)*exp(-((x-p(5))/p(6)).^2);
    %pfit = lsqcurvefit(f,p,xHist',yHist'); 

    p1 = [peak1 center1 std1];
    f = @(p,x)p(1)*exp(-((x-p(2))/(sqrt(2)*p(3))).^2);
    pfit1 = lsqcurvefit(f,p1,xHistPos',yHistPos'); 

    p2 = [peak2 center2 std2];
    pfit2 = lsqcurvefit(f,p2,xHistNeg',yHistNeg'); 

    xfit = floor(min(xHistNeg)):0.5:ceil(max(xHistPos));
    yfit1 = f(pfit1,xfit);
    line(xfit,yfit1,'Color','black','LineWidth',3,'LineStyle','-', 'Parent',axGauss);
    yfit2 = f(pfit2,xfit);
    line(xfit,yfit2,'Color','black','LineWidth',3,'LineStyle','-', 'Parent',axGauss);

    ylabel('Count','fontweight','b','fontsize',TextSize,'FontName','Palatino Linotype');
    %xlabel('Step size (nm)','fontweight','b','fontsize',12,'FontName','Palatino Linotype');
    title(GraphTitle,'fontweight','b','fontsize',TextSize,'FontName','Palatino Linotype');
    xlim([xHistNeg(NegInd)-center1/2 xHistPos(PosInd)+center1/2]);
    Text = str2mat(['Peak (-) = ' num2str(round(pfit2(2)*10)/10) ' ± ' num2str(round(pfit2(3))) ' nm'],['Peak (+) = ' num2str(round(pfit1(2)*10)/10) ' ± ' num2str(round(pfit1(3))) ' nm']);

    % Horizontal axis label (a textbox)
    %annotation('textbox',[AnnotateXLocation 0.8 0.8 0.05],'String',Text,...
    %'FontWeight','bold','FontSize',AnnotateSize,'FontName','Palatino Linotype',...
    %'FontAngle','italic','FitBoxToText','off','LineStyle','none');

    % Horizontal axis label (a textbox)
    PeakNeg = round(pfit2(2)*10)/10;
    nNeg = pfit2(1)*pfit2(3)*2.5066;
    semNeg = pfit2(3)/sqrt(nNeg);
    TextNeg = str2mat([sprintf('%1.1f',PeakNeg) ' ± ' num2str(round((semNeg*10))/10) ' nm']);
    text(PeakNeg-10,f(pfit2,PeakNeg)*1.2,TextNeg,'FontSize',TextSize,'fontweight','b')
    PeakPos = round(pfit1(2)*10)/10; 
    nPos = pfit1(1)*pfit1(3)*2.5066;
    semPos = pfit1(3)/sqrt(nPos);
    TextPos = str2mat([sprintf('%1.1f',PeakPos)  ' ± ' num2str(round((semPos*10))/10) ' nm']);
    text(PeakPos+3,f(pfit1,PeakPos)*0.98,TextPos,'FontSize',TextSize,'fontweight','b')

    print('-djpeg','-r300',[GraphTitle '.jpg']);




%f = @(x,xdata)x(1)*xdata.^2+x(2)*sin(xdata);
%x = lsqcurvefit(f,x0,xdata,ydata);
%a1*exp(-((x-b1)/c1)^2)+a2*exp(-((x-b2)/c2)^2;
    %% Modified Gaussian Fitting (Separate -ve and +ve) 
    NegShift = -25; NegVerMultiply = 1.2;
    PosShift = 5; PosVerMultiply = 1.1;
    binSpacing = 2;
    xThreshold = 0.9;
    AnnotateXLocation = 0.1;    % Range from 0 to 1. Determine horizontal location of annotation
    AnnotateSize = 15;          % Size of annotation (for center)
    TextSize = 15;
    GraphTitle = 'DDB Step Size Histogram';
    xHistPos = 0:binSpacing:ceil(max(CompiledDataPos));
    xHistNeg = floor(min(CompiledDataNeg)/binSpacing)*binSpacing:binSpacing:0;
    yHistPos=hist(CompiledDataPos,xHistPos);

    binSpace = xHistPos(2)-xHistPos(1);
    yHistNeg=hist(CompiledDataNeg,xHistNeg);

    %xHistNeg = xHist(xHist<0); yHistNeg = yHist(xHist<0);
    xHistNeg = xHistNeg(yHistNeg~=0); yHistNeg = yHistNeg(yHistNeg~=0);
    %xHistPos = xHist(xHist>=0); yHistPos = yHist(xHist>=0);
    xHistPos = xHistPos(yHistPos~=0); yHistPos = yHistPos(yHistPos~=0);
    
    % Find out variables to make the graph look full and symmetric
    CumSumNeg = cumsum(fliplr(yHistNeg)); CumSumNeg = fliplr(CumSumNeg/CumSumNeg(end));
    NegInd = find(CumSumNeg>xThreshold,1,'last');
    CumSumPos = cumsum(yHistPos); CumSumPos = CumSumPos/CumSumPos(end);
    PosInd = find(CumSumPos>xThreshold,1);
    if abs(xHistNeg(NegInd))<xHistPos(PosInd)
        while abs(xHistNeg(NegInd))<xHistPos(PosInd) && NegInd > 1
            NegInd = NegInd - 1;
        end
    elseif abs(xHistNeg(NegInd))>xHistPos(PosInd)
        while abs(xHistNeg(NegInd))>xHistPos(PosInd) && PosInd <= length(xHistPos)
            PosInd = PosInd + 1;
        end
    end

    figure;
    bar([xHistNeg(NegInd:end) xHistPos(1:PosInd)],[yHistNeg(NegInd:end) yHistPos(1:PosInd)],'r','LineWidth',1.5);
    hold on;
    %bar(xHistNeg(NegInd:end),yHistNeg(NegInd:end),'r','LineWidth',1.5);

    center1 = 8;            center2 = 10;   % center2 is for -ve peak, but the input is positive
    std1 = 10;               std2 = 20;
    peak1 = 200;             peak2 = 100;
    %p = [peak1 center1 std1 peak2 center2 std2];
    %f = @(p,x)p(1)*exp(-((x-p(2))/p(3)).^2)+p(4)*exp(-((x-p(5))/p(6)).^2);
    %pfit = lsqcurvefit(f,p,xHist',yHist'); 

    p1 = [peak1 center1 std1];
    %f = @(p,x)p(1)*exp(-((x-p(2))/(sqrt(2)*p(3))).^2);;
    f=@(p,x)(p(1)/(p(3)*sqrt(pi/2)))*sqrt(x/p(2)).*exp(-2*((x-p(2))/p(3)).^2);
    pfit1 = lsqcurvefit(f,p1,xHistPos',yHistPos'); 

    p2 = [peak2 center2 std2];
    pfit2 = lsqcurvefit(f,p2,(abs(xHistNeg))',yHistNeg'); 

    xfit1 = floor(0:0.5:ceil(max(xHistPos)));
    xfit2 = abs(floor(min(xHistNeg):0.5:0));
    yfit1 = f(pfit1,xfit1);
    plot(xfit1,yfit1,'--b','LineWidth',3);
    yfit2 = f(pfit2,xfit2);
    plot(-xfit2,yfit2,'--b','LineWidth',3);
    
    % Horizontal axis label (a textbox)
    PeakNeg = -round(pfit2(2)*10)/10;
    nNeg = integral(@(x)f(pfit2,x),0,max(abs(xHistNeg)));
    semNeg = pfit2(3)/sqrt(nNeg);
    TextNeg = str2mat([sprintf('%1.1f',PeakNeg) ' ± ' num2str(round((semNeg*10))/10) ' nm']);
    text(PeakNeg+NegShift,f(pfit2,-PeakNeg)*NegVerMultiply,TextNeg,'FontSize',TextSize,'fontweight','b')
    PeakPos = round(pfit1(2)*10)/10; 
    nPos = integral(@(x)f(pfit1,x),0,max(xHistPos));
    semPos = pfit1(3)/sqrt(nPos);
    TextPos = str2mat([sprintf('%1.1f',PeakPos)  ' ± ' num2str(round((semPos*10))/10) ' nm']);
    text(PeakPos+PosShift,f(pfit1,PeakPos)*PosVerMultiply,TextPos,'FontSize',TextSize,'fontweight','b')

    ylabel('Frequency','fontweight','b','fontsize',12,'FontName','Palatino Linotype');
    xlabel('Step size (nm)','fontweight','b','fontsize',12,'FontName','Palatino Linotype');
    title([GraphTitle ' - Modified Gaussian Fit'],'fontweight','b','fontsize',12,'FontName','Palatino Linotype');
    xlim([xHistNeg(NegInd) xHistPos(PosInd)]);

    print('-djpeg','-r300',[GraphTitle '-ModifiedGaussian.jpg']);