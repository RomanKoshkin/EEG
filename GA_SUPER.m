% This is a function that implements epoch rejection using a
% genetic algorithm. THIS IS JUST A PROOF-OF-CONCEPT implementation. We DID
% NOT use it to process the data and obtain the results reported in the 
% manuscript.

function [A, B] = GA_SUPER(filepath)
    
    % load data:
    load([filepath '/ERP.mat'])
    
    % initialize times:
    times = [-100:4:496];
    
    % GA parameters:
    TargetVar = 1;
    win0 = [find(times == -100):find(times==0)];        % window (whose mean we want to minimize)
    win = [find(times == 52):find(times==100)];         % window (whose mean we want to maximize)
    win2 = [find(times == 124):find(times==172)];       % window (whose mean we want to maximize)
    win3 = [find(times == 200):find(times==496)];       % window (whose mean we want to maximize)
    PopSize = 100;                                      % population size (too small, there's risk of premature termination
                                                        % ideally - 10% of the genePool)
    maxGens = 200;                                      % cap on the number of generations
    eliteChild = 0.1;                                   % share of children selected for mating
    nChilrenPerCouple = 10;                             % number of children per couple
    SelectionPressure = 0.9;                            % stay in the genome, others are thrown out
    
    % fitness function coefficients:
    w_P1 = 1;                                           % weight of P1
    w_N1 = 0.7;                                         % weight of N1
    w_BL = 0;                                           % weight of baseline
    w_RE = 0;                                           % weight of the rest
    
    for condition = 1:2
        clearvars -except... 
                    condition TargetVar win0 win win2 win3 PopSize...
                    maxGens eliteChild nChilrenPerCouple SelectionPressure...
                    times filepath  w_P1 w_N1 w_BL w_RE A B Cond1 Cond2...
                    generation epochNumbersInA epochNumbersInB Cnd1

        if condition == 1;
            Mat = Cond1;            % low WM load condition (as defined in ERPimg_SUPER.m)
        else
            Mat = Cond2;            % high WM load condition (as defined in ERPimg_SUPER.m)
        end
        
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


        for generation = 2:maxGens
            generation

            % eliteChild the fittest individuals (gene collections):
            fittest_ind = flipud(sortrows(fitness,1));
            fittest_ind = fittest_ind(1:end*eliteChild,:);

            % implement crossover and breeding:
            
            % Select father:
            M = fittest_ind(:,2);
            
            % we take the indices and randperm them. Select mother:
            F = randsample(fittest_ind(:,2), length(fittest_ind(:,2)));

            % a child is a vector of gene numbers (the genes are in the 
            % columns of Mat):
            count = 0;
            for i = 1:length(fittest_ind)
                for j = 1:nChilrenPerCouple % 10 children per each couple:
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

                    % define the fitness function:
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

%       uncomment if you want to plot all the epochs in the current condition:
%       for i = 1:genome_length
%           plot (times, Mat(:,i))
%           hold on
%       end
        a = 0;
        
        % create a new figure, add a relevant title:
        if ~exist('Cnd1')
            Cnd1 = figure;
            Cnd1.Name = 'Low WM load';
        else
            Cnd2 = figure;
            Cnd2.Name = 'high WM load';
        end
        
        % plot mean amplitude 200-500 ms post-stimulus:
        subplot(2,3,1)
        plot(tR)
        grid on
        title('mean amp. 200-500 poststim.')

        % plot fitness by generation:
        subplot(2,3,2)
        plot(gen_fitness)
        grid on; hold on
        title('fitness by generation')

        % plot the ERP curves before and after GA
        subplot(2,3,3)
        plot(times, mean(Mat(:,find(child(:,1))),2))
        hold on
        plot(times, mean(Mat,2), 'LineWidth', 3)
        legend ('cleaned', 'original')
        plot([times(win(1)) times(win(1))], ylim,'k-.')
        plot([times(win(end)) times(win(end))], ylim, 'k-.')
        plot([times(win0(1)) times(win0(1))], ylim,'r-.')
        plot([times(win0(end)) times(win0(end))], ylim, 'r-.')
        plot([times(win2(1)) times(win2(1))], ylim,'b-.')
        plot([times(win2(end)) times(win2(end))], ylim, 'b-.')
        str = ['Deleted Epochs=' num2str(SelectionPressure) 'TargetVar=' num2str(TargetVar)]
        Ylim=get(gca,'ylim');
        Xlim=get(gca,'xlim');
        text(Xlim(1)+5,Ylim(2)-0.2,str)
        grid on

        % show pre-stimulus variance by generation:
        subplot(2,3,4)
        plot(tvg)
        grid on; hold on
        title('pre-stim variance')
        
        % plot P1 amplitude by generation:
        subplot(2,3,5)
        plot(tP1)
        title('P1 amplitude')
        grid on; hold on
        
        % plot N1 amplitude by generation:
        subplot(2,3,6)
        plot(tN1)
        grid on; hold on
        title('N1 amplitude')
        
        % save the cleaned data to A and B to be returned by the function:
        if condition == 1;
            A = epochNumbersInA(find(child(:,1)));
        else
            B = epochNumbersInB(find(child(:,1)));
        end

    end

end