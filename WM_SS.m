% WM surround supression task

%% check parameters used
if exist([CurrDir '\Results\WMSS\' SubjID '\' SubjID '_results_sess' num2str(sess_num) '_run' num2str(run_num) '.mat'],'file')
    disp(' ');
    disp([SubjID '_run' num2str(run_num) ' has been test, please enter a new run number!']); 
    disp(' ');
    abort;      
end

results = zeros(Param.Trial.TotalNum,31);
trial_index = randperm(Param.Trial.TotalNum);
trial_index = mod(trial_index,length(Param.Stimuli.GratingOriDiff))+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% results:
% 1 trial number             %16 1st bias      
% 2 sti_loc                  %17 2nd bias      
% 3 sti_ori                  %18 trial duration        
% 4 ori of st1               %19 report1 startpoint 
% 5 ori of st2               %20 report2 startpoint
% 6 report order             %21 1st ori error
% 7 angle diff               %22 1st other ori distance   
% 8 reported ori1            %23 2nd ori error           
% 9 reported ori12           %24 2nd other ori distance
%10 trial onset              %25 WM actual diff
%11 fix onset                %26 wrong order trials
%12 ITI                      %27 1st report error
%13 st1                      %28 2nd report error
%14 ISI                      %29 delay
%15 st2                      %30 angle diff error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Go!
for trial_i = 1:Param.Trial.TotalNum

    % start
    if mod(trial_i,Param.Trial.minirun)==1  
        if trial_i == 1
            DrawFormattedText(wnd,'Press s to start!','center','center', black);
        else
            DrawFormattedText(wnd,'Take a rest! Press s to start!','center','center', black);
        end

        Screen('Flip',wnd);
        is_true = 0;
        while (is_true == 0)
            [a,b,keyCode] = KbCheck;
            if keyCode(Param.Keys.Trigger1)
                is_true = 1;
            elseif keyCode(Param.Keys.Esc)
                abort;
            end
        end
    end
    
    results(trial_i,1) = trial_i; 
    results(trial_i,2) = location_used;
    
    % define ref orientation
    temp = randperm(length(Param.Stimuli.GratingOri));
    results(trial_i,3) = Param.Stimuli.GratingOri(temp(1));

    %angle difference
    Curr_cond = trial_index(trial_i);
    %temp = randperm(lenth(Param.Stimuli.GratingOri));
    results(trial_i,7) = Param.Stimuli.GratingOriDiff(Curr_cond);

%     Curr_rotDir = ((rand-0.5)>=0)*2-1;
    % make sure the top part < 90
    Curr_rotDir = -1;
    if results(trial_i,3) > 90 
        Curr_rotDir = 1;
    end

    trial_order = ((rand-0.5)>=0);  
    if trial_order    %1ref first %0 ref second        
        results(trial_i,4) = results(trial_i,3) + (rand-0.5)*2*Param.Stimuli.OriJitter;                                        % used difference jitter
        results(trial_i,5) = results(trial_i,3) + Curr_rotDir*results(trial_i,7) + (rand-0.5)*2*Param.Stimuli.OriJitter;       % need to double check here         
    else
        results(trial_i,5) = results(trial_i,3) + (rand-0.5)*2*Param.Stimuli.OriJitter;                                        % used difference jitter
        results(trial_i,4) = results(trial_i,3) + Curr_rotDir*results(trial_i,7) + (rand-0.5)*2*Param.Stimuli.OriJitter;       % need to double check here 
    end

    results(trial_i,4) = OriCircle180(results(trial_i,4));
    results(trial_i,5) = OriCircle180(results(trial_i,5));
    
    % report order
    results(trial_i,6) = ((rand-0.5)>=0)+1;  % 1 or 2 first reported ori
    
    curr_delay = randperm(length(Param.Trial.Delay));
    results(trial_i,29) = Param.Trial.Delay(curr_delay(1));
    
    trial_onset = GetSecs;
    results(trial_i,10) = trial_onset; 

    %% ITI
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    vbl = Screen('Flip',wnd);
    results(trial_i,11) = vbl-trial_onset; 
       
    %% sti1
    [x,y]=meshgrid(-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize,-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize);
    phase = rand*2*pi; 
    angle = results(trial_i,4)/180*pi;
    sti1_final = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase));
    sti1_final(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray; %circle mask
      
    mm = Screen('MakeTexture', wnd, sti1_final);   
    Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(results(trial_i,2),1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)+Param.Stimuli.OuterSize]);
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    vbl = Screen('Flip',wnd,vbl+Param.Trial.ITI); %-Slack
    results(trial_i,12) = vbl-trial_onset-results(trial_i,11); 

    %% ISI
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    vbl = Screen('Flip',wnd,vbl+Param.Trial.StiDura); %-Slack
    results(trial_i,13) = vbl-trial_onset-sum(results(trial_i,11:12)); 

    %% sti2
    [x,y]=meshgrid(-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize,-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize);
    phase = rand*2*pi; 
    angle = results(trial_i,5)/180*pi;
    sti2_final = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase));
    sti2_final(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray; %circle mask

    mm = Screen('MakeTexture', wnd, sti2_final); 
    Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(results(trial_i,2),1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)+Param.Stimuli.OuterSize]);
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    vbl = Screen('Flip',wnd,vbl+Param.Trial.ISI); %-Slack
    results(trial_i,14) = vbl-trial_onset-sum(results(trial_i,11:13));

    %% delay
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    vbl = Screen('Flip',wnd,vbl+Param.Trial.StiDura);     %-Slack   
    results(trial_i,15) = vbl-trial_onset-sum(results(trial_i,11:14));

    %% report 1
    [x,y]=meshgrid(-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize,-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize);
    phase = rand*2*pi; 
    angle = rand*pi;
    results(trial_i,19) = angle*180/pi;

    sti3_final = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase));
    sti3_final(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray; %circle mask
      
    mm = Screen('MakeTexture', wnd, sti3_final);   
    Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(results(trial_i,2),1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)+Param.Stimuli.OuterSize]);
    
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    if results(trial_i,6) == 1
        DrawFormattedText(wnd,'1','center',Param.Settings.ScrnResolution(4)/2-Param.Fixation.Offset,Param.Fixation.FontColor);
    else
        DrawFormattedText(wnd,'2','center',Param.Settings.ScrnResolution(4)/2-Param.Fixation.Offset,Param.Fixation.FontColor);
    end
    
%     Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    vbl = Screen('Flip',wnd,vbl + results(trial_i,29));     %-Slack 

    is_true = 0;
    Angle_temp = 0;
    while (is_true == 0)
        [a,b,keyCode] = KbCheck; 
        if keyCode(Param.Keys.Right)
            Angle_temp = Angle_temp+Param.Trial.ResponsePrec; 
        elseif keyCode(Param.Keys.Left)
            Angle_temp = Angle_temp-Param.Trial.ResponsePrec;
        elseif keyCode(Param.Keys.Space)
            is_true = 1;
        elseif keyCode(Param.Keys.Esc)
            abort;
        end
        
        Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(results(trial_i,2),1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)+Param.Stimuli.OuterSize],Angle_temp);    
        Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
        if results(trial_i,6) == 1
            DrawFormattedText(wnd,'1','center',Param.Settings.ScrnResolution(4)/2-Param.Fixation.Offset, Param.Fixation.FontColor);
        else
            DrawFormattedText(wnd,'2','center',Param.Settings.ScrnResolution(4)/2-Param.Fixation.Offset, Param.Fixation.FontColor);
        end        
%         Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
        Screen('Flip',wnd);
    end

    results(trial_i,8) = OriCircle180(Angle_temp+angle/pi*180);

    %% ISI
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    vbl = Screen('Flip',wnd); %-Slack
       
    %% report 2
    [x,y]=meshgrid(-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize,-Param.Stimuli.OuterSize:Param.Stimuli.OuterSize);
    phase = rand*2*pi; 
    angle = rand*pi;
    results(trial_i,20) = angle*180/pi;

    sti4_final = gray*(1+Param.Stimuli.GratingContrast*cos(2*pi*Param.Stimuli.Spatial_freq*(y*sin(angle)+x*cos(angle))+phase));
    sti4_final(sqrt(x.^2 + y.^2) > Param.Stimuli.OuterSize) = gray; %circle mask
      
    mm = Screen('MakeTexture', wnd, sti4_final);   
    Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(results(trial_i,2),1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)+Param.Stimuli.OuterSize]);
 
    Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
    if results(trial_i,6) == 1
        DrawFormattedText(wnd,'2','center',Param.Settings.ScrnResolution(4)/2-Param.Fixation.Offset, Param.Fixation.FontColor);
    else
        DrawFormattedText(wnd,'1','center',Param.Settings.ScrnResolution(4)/2-Param.Fixation.Offset, Param.Fixation.FontColor);
    end
%     Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
    vbl = Screen('Flip',wnd,vbl+Param.Trial.ISI); %-Slack

    is_true = 0;
    Angle_temp = 0;
    while (is_true == 0)
        [a,b,keyCode] = KbCheck;
        if keyCode(Param.Keys.Right)
            Angle_temp = Angle_temp+Param.Trial.ResponsePrec;     
        elseif keyCode(Param.Keys.Left)
            Angle_temp = Angle_temp-Param.Trial.ResponsePrec;
        elseif keyCode(Param.Keys.Space)
            is_true = 1;
        elseif keyCode(Param.Keys.Esc)
            abort;
        end

        Screen('DrawTexture', wnd, mm,[],[Param.Stimuli.Locations(results(trial_i,2),1)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)-Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),1)+Param.Stimuli.OuterSize,Param.Stimuli.Locations(results(trial_i,2),2)+Param.Stimuli.OuterSize],Angle_temp);    
        Screen(wnd,'FillOval', black, CenterRect([0 0 Param.Fixation.OvalSize Param.Fixation.OvalSize],Param.Settings.ScrnResolution));
        if results(trial_i,6) == 1
            DrawFormattedText(wnd,'2','center',Param.Settings.ScrnResolution(4)/2-Param.Fixation.Offset, Param.Fixation.FontColor);
        else
            DrawFormattedText(wnd,'1','center',Param.Settings.ScrnResolution(4)/2-Param.Fixation.Offset, Param.Fixation.FontColor);
        end        
%         Screen('DrawLines', wnd, Param.Fixation.CrossLoc, Param.Fixation.CrossWidth, Param.Fixation.CrossColor, [], 1);  
        Screen('Flip',wnd);
    end

    results(trial_i,9) = OriCircle180(Angle_temp+angle/pi*180);  

    % calculate bias
    if results(trial_i,6) == 1  
        results(trial_i,21) = CalAngleDiff90(abs(results(trial_i,4)-results(trial_i,8)));       
        results(trial_i,22) = CalAngleDiff90(abs(results(trial_i,5)-results(trial_i,8)));
        results(trial_i,23) = CalAngleDiff90(abs(results(trial_i,5)-results(trial_i,9)));
        results(trial_i,24) = CalAngleDiff90(abs(results(trial_i,4)-results(trial_i,9)));
    else
        results(trial_i,21) = CalAngleDiff90(abs(results(trial_i,4)-results(trial_i,9)));       
        results(trial_i,22) = CalAngleDiff90(abs(results(trial_i,5)-results(trial_i,9)));
        results(trial_i,23) = CalAngleDiff90(abs(results(trial_i,5)-results(trial_i,8)));
        results(trial_i,24) = CalAngleDiff90(abs(results(trial_i,4)-results(trial_i,8)));
    end

    results(trial_i,25) = CalAngleDiff90(abs(results(trial_i,4)-results(trial_i,5)));
    results(trial_i,30) = CalAngleDiff90(abs(results(trial_i,8)-results(trial_i,9)));
    results(trial_i,31) = results(trial_i,30)-results(trial_i,25);

    if (results(trial_i,21) + results(trial_i,22)) - results(trial_i,25) > 0.000001
        results(trial_i,16) = results(trial_i,21);
    else
        results(trial_i,16) = -results(trial_i,21);
    end

    if (results(trial_i,23) + results(trial_i,24)) - results(trial_i,25) > 0.000001
        results(trial_i,17) = results(trial_i,23);
    else
        results(trial_i,17) = -results(trial_i,23);
    end

    % remember wrong order
    if (results(trial_i,22) < results(trial_i,21)) && (results(trial_i,24) < results(trial_i,23))
        results(trial_i,26) = 1;
    end

    % remember wrong order
    if results(trial_i,8)>=results(trial_i,4) && results(trial_i,8)>=results(trial_i,5) && results(trial_i,9)>=results(trial_i,5) && results(trial_i,9)>=results(trial_i,4)
        results(trial_i,26) = 1;
    end

    % remember wrong order
    if results(trial_i,8)<=results(trial_i,4) && results(trial_i,8)<=results(trial_i,5) && results(trial_i,9)<=results(trial_i,5) && results(trial_i,9)<=results(trial_i,4)
        results(trial_i,26) = 1;
    end

    if results(trial_i,6) == 1 
        results(trial_i,27) = results(trial_i,16);
        results(trial_i,28) = results(trial_i,17);
    else
        results(trial_i,27) = results(trial_i,17);
        results(trial_i,28) = results(trial_i,16);
    end

    results(trial_i,18) = GetSecs - results(trial_i,10);
end
       
%% save the data
resultsDir = [CurrDir '\Results\WMSS\' SubjID '\'];
if ~isdir(resultsDir)                                                       %isder和mkder可以用来创建新文件夹
    mkdir(resultsDir);
end

cd(resultsDir);
results_name = [SubjID '_results_sess' num2str(sess_num) '_run' num2str(run_num) '.mat'];
save(results_name,'results','Param');
cd(CurrDir); 


%% plot
% distribution
repulsion_mean = zeros(4,length(Param.Stimuli.GratingOriDiff));
repulsion_std  = zeros(4,length(Param.Stimuli.GratingOriDiff));
figure(1);
for condi = 1:length(Param.Stimuli.GratingOriDiff)
    curr_condi = Param.Stimuli.GratingOriDiff(condi);
    repulsion_mean(1,condi) = mean(results((results(:,7)==curr_condi & (~results(:,26)) & (results(:,29)==Param.Trial.Delay(1))),16));
    repulsion_mean(2,condi) = mean(results((results(:,7)==curr_condi & (~results(:,26)) & (results(:,29)==Param.Trial.Delay(1))),17));
    repulsion_std(1,condi)  = std (results((results(:,7)==curr_condi & (~results(:,26)) & (results(:,29)==Param.Trial.Delay(1))),16));
    repulsion_std(2,condi)  = std (results((results(:,7)==curr_condi & (~results(:,26)) & (results(:,29)==Param.Trial.Delay(1))),17));

    repulsion_mean(3,condi) = mean(results((results(:,7)==curr_condi & (~results(:,26)) & (results(:,29)==Param.Trial.Delay(1))),27));
    repulsion_mean(4,condi) = mean(results((results(:,7)==curr_condi & (~results(:,26)) & (results(:,29)==Param.Trial.Delay(1))),28));
    repulsion_std(3,condi)  = std (results((results(:,7)==curr_condi & (~results(:,26)) & (results(:,29)==Param.Trial.Delay(1))),27));
    repulsion_std(4,condi)  = std (results((results(:,7)==curr_condi & (~results(:,26)) & (results(:,29)==Param.Trial.Delay(1))),28));
    
    subplot(4,length(Param.Stimuli.GratingOriDiff),condi+length(Param.Stimuli.GratingOriDiff)*0);hist(results((results(:,7)==curr_condi & (~results(:,26)) & (results(:,29)==Param.Trial.Delay(1))),16));title(['ori1-cond' num2str(condi)]);
    subplot(4,length(Param.Stimuli.GratingOriDiff),condi+length(Param.Stimuli.GratingOriDiff)*1);hist(results((results(:,7)==curr_condi & (~results(:,26)) & (results(:,29)==Param.Trial.Delay(1))),17));title(['ori2-cond' num2str(condi)]);
    subplot(4,length(Param.Stimuli.GratingOriDiff),condi+length(Param.Stimuli.GratingOriDiff)*2);hist(results((results(:,7)==curr_condi & (~results(:,26)) & (results(:,29)==Param.Trial.Delay(1))),27));title(['rep1-cond' num2str(condi)]);
    subplot(4,length(Param.Stimuli.GratingOriDiff),condi+length(Param.Stimuli.GratingOriDiff)*3);hist(results((results(:,7)==curr_condi & (~results(:,26)) & (results(:,29)==Param.Trial.Delay(1))),28));title(['rep2-cond' num2str(condi)]);
end

% error
figure(2);
subplot(1,2,1);
errorbar(Param.Stimuli.GratingOriDiff,repulsion_mean(1,1:length(Param.Stimuli.GratingOriDiff)),repulsion_std(1,1:length(Param.Stimuli.GratingOriDiff)),'-.ob','MarkerSize',10);
hold on;
errorbar(Param.Stimuli.GratingOriDiff,repulsion_mean(2,1:length(Param.Stimuli.GratingOriDiff)),repulsion_std(2,1:length(Param.Stimuli.GratingOriDiff)),'-.*r','MarkerSize',10);
axis([-10 60 -30 30]);
title('delay = 3s; ori');
xlabel('Orientation Difference', 'FontSize',14);
ylabel('Attraction vs. Repulsion', 'FontSize',14);
legend('1st Orientaion','2nd Orientaion','FontSize',14);

subplot(1,2,2);
errorbar(Param.Stimuli.GratingOriDiff,repulsion_mean(3,1:length(Param.Stimuli.GratingOriDiff)),repulsion_std(3,1:length(Param.Stimuli.GratingOriDiff)),'-.ob','MarkerSize',10);
hold on;
errorbar(Param.Stimuli.GratingOriDiff,repulsion_mean(4,1:length(Param.Stimuli.GratingOriDiff)),repulsion_std(4,1:length(Param.Stimuli.GratingOriDiff)),'-.*r','MarkerSize',10);
axis([-10 60 -30 30]);
title('delay = 3s; report');
xlabel('Orientation Difference', 'FontSize',14);
ylabel('Attraction vs. Repulsion', 'FontSize',14);
legend('1st Report','2nd Report','FontSize',14);

% angle diff error
figure(3);
repulsion2_mean = zeros(4,length(Param.Stimuli.GratingOriDiff));
repulsion2_std  = zeros(4,length(Param.Stimuli.GratingOriDiff));
for condi = 1:length(Param.Stimuli.GratingOriDiff)
    curr_condi = Param.Stimuli.GratingOriDiff(condi);
    repulsion2_mean(1,condi) = mean(results((results(:,7)==curr_condi & (results(:,29)==Param.Trial.Delay(1))),31));  
    repulsion2_std(1,condi)  = std (results((results(:,7)==curr_condi & (results(:,29)==Param.Trial.Delay(1))),31));    
end

errorbar(Param.Stimuli.GratingOriDiff,repulsion2_mean(1,1:length(Param.Stimuli.GratingOriDiff)),repulsion2_std(1,1:length(Param.Stimuli.GratingOriDiff)),'-.ob','MarkerSize',10);
axis([-10 60 -30 30]);
title('delay = 3s; report');
xlabel('Orientation Difference', 'FontSize',14);
ylabel('Attraction vs. Repulsion', 'FontSize',14);

% angle diff error
figure(4);
scatter(results(:,25),results(:,31));
axis([-10 60 -30 30]);
title('reported error vs actual ori');

figure(5);
hist(results(:,16));
title('error distribution of st1');

figure(6);
hist(results(:,17));
title('error distribution of st2');
%%
reset_test_gamma;
ShowCursor;
Screen('CloseAll');
 
delete *.asv
