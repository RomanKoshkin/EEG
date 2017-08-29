function [A, B] = GA_SUPER(filepath)
for condition = 1:2
clearvars -except condition filepath A B
% weights:
w_P1 = 1;
w_N1 = 0.7;
w_BL = 0;
w_RE = 0;

load([filepath '/ERP.mat'])
times = [-100:4:396];

if condition == 1;
    Mat = Cond1;
else
    Mat = Cond2;
end
%%%%%% PARAMETERS:%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TargetVar = 1;
win0 = [find(times == -100):find(times==0)];
win = [find(times == 52):find(times==100)];         % winow (whose mean we want to maximize)
win2 = [find(times == 124):find(times==172)];       % winow (whose mean we want to maximize)
win3 = [find(times == 200):find(times==396)];       % winow (whose mean we want to maximize)

PopSize = 100;       % population size (too small, there's risk of premature termination
                     % ideally - 10% of the genePool)
nGens = 300;         % number of generations
eliteChild = 0.1;    % share of children selected for mating
nChilrenPerCouple = 10;
SelectionPressure = 0.9; % stay in the genome, others are thrown out
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

genePool = size(Mat,2);      % gene pool (from which to select)
genome_length = round(genePool*SelectionPressure); % number of genes that will stay in the clean species     

fitness = zeros(PopSize,1);

% generate the initial population such that each individual has only 80%
% of the genes randomly sampled from the gene pool:
tic
for ind = 1:PopSize
    temp = [1:genePool]';
    temp(randsample(genePool, genePool-genome_length)) = 0; %we zero out 20% of the gene nubers randomly
    population(:,ind) = temp; % create a population based on 80% of the gene pool
    
    tempP1 =  mean(mean(Mat(win,  find(temp>0)), 2));
    tempN1 =  mean(mean(Mat(win2, find(temp>0)), 2));
    tempVar = var (mean(Mat(win0, find(temp>0)), 2));
    tempRest =  mean(mean(Mat(win3, find(temp>0)), 2));

    fitness(ind,1:2) = [(w_P1*tempP1-w_N1*tempN1) ind];
    %fitness(ind,1:2) = [(w_P1*tempP1-w_N1*tempN1)/(w_RE*abs(tempRest)+w_BL*tempVar) ind];
end
gen_fitness(1) = mean(fitness(:,1)); %compute mean fitness of the generation
tvg(1) = tempVar;
tP1(1) = tempP1;
tN1(1) = tempN1;
tR(1)  = tempRest;


for generation = 2:nGens
generation

% eliteChild the fittest individuals (gene collections):
fittest_ind = flipud(sortrows(fitness,1));
fittest_ind = fittest_ind(1:end*eliteChild,:);

% implement crossover and breeding:
M = fittest_ind(:,2);                                   % select father
% we take the indices and randperm them
F = randsample(fittest_ind(:,2), length(fittest_ind(:,2))); % select mother

% child is a vector with gene numbers (you get genes from the col of Mat:
count = 0;
for i = 1:length(fittest_ind)
    for j = 1:nChilrenPerCouple % TEN children per each couple:
        count = count + 1;
        temp = [1:genePool]'; temp(randsample(genePool, genePool-genome_length)) = 0;
        xover_mat = temp>0; % random parent gene mixing matrix
        inv_xover_mat = ~temp>0;
       
        child(:,count) = double(or(xover_mat.*population(:,M(i)), inv_xover_mat.*population(:, F(i))));
        
        while sum(child(:,count)) > genome_length
            child(randsample(find(child(:,count)==1),1), count) = 0;
        end
        
        while sum(child(:,count)) < genome_length
            child(randsample(find(child(:,count)==0),1), count) = 1;
        end
        
% FITNESS FUNCTION:
        tempP1 = mean(mean(Mat(win, find(child(:,count))),2));
        tempN1 = mean(mean(Mat(win2,find(child(:,count))),2));
        tempVar = var(mean(Mat(win0,find(child(:,count))),2));
        tempRest = mean(mean(Mat(win3,find(child(:,count))),2));
       
                tvg(generation) = tempVar; % for plotting later
                tP1(generation) = tempP1;
                tN1(generation) = tempN1;
                tR(generation) = tempRest;
        child_fitness(count,1:2) = [(w_P1*tempP1-w_N1*tempN1) count];
        %child_fitness(count,1:2) = [(w_P1*tempP1-w_N1*tempN1)/(w_RE*abs(tempRest)+w_BL*tempVar) count];
% END FITNESS FUNCTION
    end
end
gen_fitness(generation) = mean(child_fitness(:,1));
if gen_fitness(generation) == gen_fitness(generation - 1)
    break
end
fitness = child_fitness;
population = child;
end
toc

% subplot(2,3,1)
% plot(tR)
% grid on
% title('mean amp. 200-400 poststim.')

% % for i = 1:genome_length
% %     plot (times, Mat(:,i))
% %     hold on
% % end
a = 0;


% % 
% % plot(mean(Mat(:,:),2), 'LineWidth', 2)
% % plot(mean(Mat(:,find(child(:,1))),2), 'LineWidth', 3)
% % 
% % plot([times(win(1)) times(win(1))], Ylim,'k-.')
% % plot([times(win(end)) times(win(end))], Ylim, 'k-.')
% % plot([times(win0(1)) times(win0(1))], Ylim,'r-.')
% % plot([times(win0(end)) times(win0(end))], Ylim, 'r-.')
% % plot([times(win2(1)) times(win2(1))], Ylim,'b-.')
% % plot([times(win2(end)) times(win2(end))], Ylim, 'b-.')

% subplot(2,3,2)
% plot(gen_fitness)
% grid on; hold on
% title('fitness')
% 
% subplot(2,3,3)
% plot(times, mean(Mat(:,find(child(:,1))),2))
%     
% hold on
% plot(times, mean(Mat,2), 'LineWidth', 3)
% legend ('cleaned', 'original')
% plot([times(win(1)) times(win(1))], ylim,'k-.')
% plot([times(win(end)) times(win(end))], ylim, 'k-.')
% plot([times(win0(1)) times(win0(1))], ylim,'r-.')
% plot([times(win0(end)) times(win0(end))], ylim, 'r-.')
% plot([times(win2(1)) times(win2(1))], ylim,'b-.')
% plot([times(win2(end)) times(win2(end))], ylim, 'b-.')
% str = ['Deleted Epochs=' num2str(SelectionPressure) 'TargetVar=' num2str(TargetVar)]
% Ylim=get(gca,'ylim');
% Xlim=get(gca,'xlim');
% text(Xlim(1)+5,Ylim(2)-0.2,str)
% grid on
% 
% subplot(2,3,4)
% plot(tvg)
% grid on; hold on
% title('pre-stim variance')
% subplot(2,3,5)
% plot(tP1)
% title('P1 amplitude')
% grid on; hold on
% subplot(2,3,6)
% plot(tN1)
% grid on; hold on
% title('N1 amplitude')

if condition == 1;
    A = epochNumbersInA(find(child(:,1)));
else
    B = epochNumbersInB(find(child(:,1)));
end

end

end