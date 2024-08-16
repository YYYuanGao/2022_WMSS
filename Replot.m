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