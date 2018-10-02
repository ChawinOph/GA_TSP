classdef SimpleGA < handle
    %SimpleGA Summary of this class goes here
    %   This class is designed for a TSP
    
    properties (Constant)
        n_bit = 1000; % chormosome length
    end
    
    properties
        filename            % file name of the given set of points
        n_pop               % number of population
        pool                % n_bit x n_pop array of current chormosomes
        p_c                 % cross over prob (single crossover)
        p_m                 % mutation prob
        n_crossover         % number of crossover indice of parent chromosomes (no. of cutting lines)
        n_mutation          % number of mutation indice in each parent chromosome
        n_iter              % number of evaluations/n_pop
        trunc_rate          % truncation selection ratio to improve diversity (1/n fraction form only!)
        b_find_longest      % boolean: choose goal
        n_tournament        % number of parents for competing before geeting only top two
    end
    
    properties % Dependent variables
        fitness         % 1 x n_pop:  array of fitness correspoinding to each chromosome
        points          % 1000 x 2 set of points from the text file
        init_points     % fixed 1000 x 2 set of points from the text file
        x_mat           % n_bit x n_pop array: repeated n_pop columns of x position of points
        y_mat           % n_bit x n_pop array: repeated n_pop columns of y position of points
        fittest         % 1 x n_eval:  array of bbest fitness value of each iteration
        x_mat_trunc     % n_bit x (n_pop/trunc_rate) array: repeated (n_pop/trunc_rate) columns of x position of points
        y_mat_trunc     % n_bit x (n_pop/trunc_rate) array: repeated (n_pop/trunc_rate) columns of y position of points
        fitness_hist    % n_iter x n_pop:  array of fitness correspoinding to each chromosome over all iterations
    end
    
    methods
        %% Constructor
        function this = SimpleGA(filename, n_pop, p_c, p_m, n_crossover, n_mutation, n_eval, n_tournament, trunc_rate, b_find_longest)
            %SimpleGA: Construct an instance of this class
            %   Detailed explanation goes here
            this.filename = filename;
            this.n_pop = n_pop;  this.p_c = p_c;  this.p_m = p_m;
            this.n_crossover = n_crossover; this.n_mutation = n_mutation;
            this.points = this.importPath(this.filename);
            this.init_points = this.importPath(this.filename);
            this.x_mat = repmat(this.points(:,1),1,this.n_pop);
            this.y_mat = repmat(this.points(:,2),1,this.n_pop);
            this.trunc_rate = trunc_rate;
            this.x_mat_trunc = repmat(this.points(:,1),1,this.n_pop/this.trunc_rate);
            this.y_mat_trunc = repmat(this.points(:,2),1,this.n_pop/this.trunc_rate);
            this.fittest = zeros(n_eval, 1); % the same size as total eval
            this.fitness_hist = zeros(n_eval, n_pop); % the same size as total eval
            this.n_iter = n_eval/n_pop; % according to the populaiton size
            this.b_find_longest = b_find_longest;
            this.n_tournament = n_tournament;
            
                      
            % Create a pool of chromosomes
            this.pool = rand(this.n_bit, this.n_pop);
            this.updateFitness()
        end
        
        %% Getters
        function fitness = getFitness(this)
            fitness = this.fitness;
        end
        
        function fittest = getFittest(this)
            fittest = this.fittest;
        end
        
        %% Member Functions
        function evaluate(this)
            % evaluate: Perform an evaluation (evolving to get a new set of
            % n_pop offspring)
            for n = 1 : this.n_iter
                % create a new array of offspring
                offspring = zeros(this.n_bit, 2*ceil(this.n_pop/2)*(1/this.trunc_rate - 1));
                
                for n_t = 1 : (1/this.trunc_rate - 1)
                    for i = 1 : ceil(this.n_pop/2)                      
                        % Tournament Selection
                        % Select a pair of parent chromosome (simultaneously) base on the fitness
                        if this.b_find_longest % choose by weighting the lengths
                            candidate_parent_indice = randsample(1:this.n_pop, this.n_tournament, true, this.fitness);
                            [~, sorted_candidate_indc] = sort(this.fitness(candidate_parent_indice), 'descend');
                        else % choose by weighting reciprocal of the lengths
                            candidate_parent_indice = randsample(1:this.n_pop, this.n_tournament, true, 1./this.fitness);
                            [~, sorted_candidate_indc] = sort(1./this.fitness(candidate_parent_indice), 'descend');
                        end
                        parent_indice = candidate_parent_indice(sorted_candidate_indc(1:2));
                        parents = this.pool(:, parent_indice);
                        
                        % CROSSOVER
                        if rand <= this.p_c
                            % Crossover the parent at a randomly chosen point
                            % The index starting from 1 to crossover_indx is the
                            % portion that will be kept in the first parent
                            if this.n_crossover == 1
                                cross_start_indx = randperm(this.n_bit - 1, 1) + 1; % - 1 to prevent from selecting the last bit
                                swapped_bits = parents(cross_start_indx : end, 2);
                                parents(cross_start_indx : end, 2) = parents(cross_start_indx : end, 1);
                                parents(cross_start_indx : end, 1) = swapped_bits;
                            else
                                cross_indc = randperm(this.n_bit - 1, this.n_crossover);% -1 to prevent from selecting the last bit
                                cross_indc = sort(cross_indc, 'ascend');
                                for j = 1:length(cross_indc)
                                    % crossover only after the odd crossover
                                    % point until reaching the next even crossover point
                                    if rem(j,2) == 1
                                        cross_start_indx = cross_indc(j) + 1;
                                        if j == length(cross_indc)
                                            % end_index is the last bit for the last
                                            % crossover after the odd point
                                            cross_end_indx = this.n_bit;
                                        else
                                            cross_end_indx = cross_indc(j + 1);
                                        end
                                        swapped_bits = parents(cross_start_indx : cross_end_indx, 2);
                                        parents(cross_start_indx : cross_end_indx, 2) = parents(cross_start_indx : cross_end_indx, 1);
                                        parents(cross_start_indx : cross_end_indx, 1) = swapped_bits;
                                    end
                                end
                            end
                        end
                        
                        % MUTATION
                        if rand <= this.p_m
                            % Choose an index for mutation
                            mutation_indc = randperm(this.n_bit, this.n_mutation);
                            for mutation_indx = mutation_indc
                                parents(mutation_indx, 1) = 1 - parents(mutation_indx, 1);
                                parents(mutation_indx, 2) = 1 - parents(mutation_indx, 2);
                            end
                        end
                        
                        % Copy the parents as two new offspring to the next gen
                        % population group
                        offspring(:, 2*i-1) = parents(:, 1);
                        offspring(:, 2*i) = parents(:, 2);
                    end
                end
                
                if rem(this.n_pop, 2) == 1
                    % If n is odd, one new population member can be discarded at random
                    discard_indx =  randperm(length(offspring),(1/this.trunc_rate - 1));
                    offspring = offspring(1:end ~= discard_indx);
                end
                
                % Update the pool with a new set of offspring and parent by getting
                % only the top 100*trunc_rate% of the combined pool
                combined_pool = [this.pool, offspring];
                combined_fitness = this.calcFitness(combined_pool);
                
                if this.b_find_longest % keep longer paths
                    [~, sorted_fitness_order] = sort(combined_fitness,'descend');
                else % keep shorter paths
                    [~, sorted_fitness_order] = sort(combined_fitness,'ascend');
                end
                
                combined_pool = combined_pool(:,sorted_fitness_order);
                % Get the 100*trunc_rate%
                this.pool = combined_pool(:, 1:this.n_pop);
                % Update the fitness again after the truncation
                this.fitness = combined_fitness(:, 1:this.n_pop);
                
                % Record data as well as getting the fittest index (default is 1, but just in case)
                if this.b_find_longest
                    [current_fittest, fittest_indx] = max(this.fitness);
                else
                    [current_fittest, fittest_indx] = min(this.fitness);
                end
                % Store fitness values of every chromosome in the current
                % iteration
                start_inx = 1 + (n - 1)*this.n_pop;
                this.fittest(start_inx : start_inx + this.n_pop - 1) = current_fittest*ones(this.n_pop, 1);
                this.fitness_hist(start_inx : start_inx + this.n_pop - 1 , :) = repmat(this.fitness, this.n_pop, 1);            
                [ ~, fittest_path_order] = sort(this.pool(:, fittest_indx), 'descend');
                
                % Sort the rows of the given set of points based on the fittest value 
                % and sort all choromosome in the same way
                % This way will give more chance of getting better genes in
                % bad chromosome and tightening the linkage because bits
                % that are supposed to be near each other will have less
                % chance to break down
                this.x_mat = this.x_mat(fittest_path_order, :);
                this.y_mat = this.y_mat(fittest_path_order, :);
                this.x_mat_trunc = this.x_mat_trunc(fittest_path_order, :);
                this.y_mat_trunc = this.y_mat_trunc(fittest_path_order, :);
                this.points = this.points(fittest_path_order, :);
                this.pool = this.pool(fittest_path_order, :);            
            end
        end
        
        function fitness = calcFitness(this, offspring)
            % calcFitness: Return 1 x n_offspring array of fitness values
            [~, sorted_order] = sort(offspring,'descend');
            
            % sort each column of the x_mat or y_mat based on the
            % each corresponding column of the sorted_order matrix
            [m, n]= size(sorted_order);
            sorted_x_mat = this.x_mat_trunc(sub2ind([m n],sorted_order,repmat(1:n,m,1)));
            sorted_y_mat = this.y_mat_trunc(sub2ind([m n],sorted_order,repmat(1:n,m,1)));
            
             % Append the first point at the bottom
            traveled_x_mat = [sorted_x_mat; sorted_x_mat(1,:)];
            traveled_y_mat = [sorted_y_mat; sorted_y_mat(1,:)] ;
            
            % Calculate path length from each column of the sorted order
            % matrix
            fitness = sum(sqrt(diff(traveled_x_mat).^2 + diff(traveled_y_mat).^2));     
        end
        
        function updateFitness(this)
            % updateFitness: Update fitness array
            [~, sorted_order] = sort(this.pool,'descend');
            
            % sort each column of the x_mat or y_mat based on the
            % each corresponding column of the sorted_order matrix
            [m, n]= size(sorted_order);
            sorted_x_mat = this.x_mat(sub2ind([m n],sorted_order,repmat(1:n,m,1)));
            sorted_y_mat = this.y_mat(sub2ind([m n],sorted_order,repmat(1:n,m,1)));
            
            % Append the first point at the bottom
            traveled_x_mat = [sorted_x_mat; sorted_x_mat(1,:)];
            traveled_y_mat = [sorted_y_mat; sorted_y_mat(1,:)];
            
            % Calculate path length from each column of the sorted order
            % matrix
            this.fitness = sum(sqrt(diff(traveled_x_mat).^2 + diff(traveled_y_mat).^2));
        end
        
        function points = importPath(~, text_filename)
            % Import the text file or array
            if ischar(text_filename)
                dlist = dir(text_filename);
                file_directory = [dlist.folder '\' dlist.name] ;
                points = dlmread(file_directory);
            elseif isfloat(text_filename)
                points = text_filename;
            end
        end
        
        function plotSemiLog(this, color)
            plot(1 : this.n_iter*this.n_pop, this.getFittest(), color)
        end
         
    end
end

