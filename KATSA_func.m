function [gbest,gbestval,convergence,time] = KATSA_func(fhd,D,Particle_Number,Max_Gen,dmin,dmax,varargin)
    tic; % Start stopwatch timer
    N=Particle_Number;
    seeds_lb=ceil(N*0.1);% upper limit of seed's number
    seeds_ub=ceil(N*0.25);% lower limit of seed's number
    best_params=zeros(1,D);
    best=inf;
    convergence=zeros(1,Max_Gen);
    ST = 0.1*ones(1,N);
    ctx = zeros(1,N);
    % JJH: 设置 distanceMatrix
    treeDistance = zeros(N,N);
    nk = 2; %最近的5棵树
    flin= zeros(1,N);
    jinlin = 0;
    fjinlin = 0;
    qianyi = 0;
    fqianyi = 0;
    B_neighbour_set = zeros(nk+1,D);
    trees=zeros(N,D);    % tree population
    obj_trees=zeros(1,N);      % objective function value for each tree
    tic;
    %% Initialization
    for i=1:N
        for j=1:D
            Boundary_no= size(dmax,2);
            if Boundary_no==1
                trees(i,j)=dmin+(dmax-dmin)*rand;
            end
            if Boundary_no>1
                for k=1:D
                    dmax_k=dmax(k);
                    dmin_k=dmin(k);
                    trees(:,k)=rand(N,1).*(dmax_k-dmin_k)+dmin_k;
                end
            end
        end
    end
    obj_trees=feval(fhd,trees',varargin{:});
    BNN_tree.distance = zeros(nk,D);
    BNN_tree.index = zeros(nk,1);
    fes=N;              % is the fes counter.
    for iter=1:Max_Gen
        for i=1:N
            Number_of_seeds=fix(seeds_lb+(seeds_ub-seeds_lb)*rand)+1;
            seeds=zeros(Number_of_seeds,D);
            obj_seeds=zeros(1,Number_of_seeds);
            %Find the best one of trees
            [best,min_indis]=min(obj_trees);
            %% set k neighbour tree
            best_params=trees(min_indis,:);%best tree
            treeDistance = pdist2(best_params,trees); %2. 计算距离
            [A, B] = sort(treeDistance);
            BNN_tree.distance = A(2:nk+1);
            BNN_tree.index =B(2:nk+1);
            tBNN_tree.distance = A(nk+2:N);%不是最近的5个（其他24个）
            tBNN_tree.index =B(nk+2:N);
            B_neighbour_set(1,:) = best_params;
            ss = 2;
            for x = BNN_tree.index
                B_neighbour_set(ss,:) = trees(x,:);
                ss = ss + 1;
            end
            %% whether the tree in the range of best_params' neighbour
            if all(ismember(trees(i,:),B_neighbour_set))%如果这棵树在邻近区域
                jinlin = jinlin + 1;
                flin(i) = 0;
                ST(i) = 0.8;
            else%如果树不在best_params附近
                flin(i) = flin(i) + 1;
                fjinlin = fjinlin + 1;
                if(flin(i) >= 10)%如果连续10次搜索都不在附近 那么迁移树
                    qianyi = qianyi + 1;
                    tnk=fix(rand*(N-nk-1))+1;
                    tn=tBNN_tree.index(tnk);%选树(1~24)
                    rnk=fix(rand*nk)+1;%近的树2
                    r=BNN_tree.index(rnk);
                    for d = 1:D
                        trees(i,d) = best_params(d)+(trees(tn,d)-trees(r,d))*cos(rand*2*pi);%Eq.8
                    end
                    Flag4ub=trees(i,:)>dmax;
                    Flag4lb=trees(i,:)<dmin;
                    trees(i,:)=(trees(i,:).*(~(Flag4ub+Flag4lb)))+dmax.*Flag4ub+dmin.*Flag4lb;
                    obj_trees(i)=feval(fhd,trees(i,:)',varargin{:});
                    ST(i) = 0.8;%迁移后在局部搜索
                    flin(i) = 0;
                else%没有10次 继续搜索
                    fqianyi= fqianyi + 1;
                    ST(i) = 0.2;%没迁移就全局搜索
                end%if(flin(i) >= 10)
            end%if(all(ismember(trees(i,:),B_neighbour_set)))%如果连续10次搜索都不在附近 那么迁移树 并进行一次局部搜索
            %% seed iteration 
            if ctx(i) >= 1%如果这棵树在上一轮发生了替换 继续搜索邻近
                % line:124
                ctx(i) = 0;
                ST(i) = 0.8;%3点改
                for j=1:Number_of_seeds
                    r=fix(rand*N)+1;
                    while(i==r)
                        r=fix(rand*N)+1;
                    end
                    for d=1:D
                        if(rand<ST(i))
                            seeds(j,d)=trees(i,d)+(trees(i,d)-trees(r,d))*sin(rand*pi*2)*0.05;% Eq.10
                            %种子应该在树的非常附近播撒
                        else
                            seeds(j,d)=(trees(i,d)+best_params(d))/2+(trees(tn,d)-best_params(d))*cos(rand*pi*2); %Eq.12
                            %如果是一个欺骗位置 依靠最佳树跳出
                        end
                    end;
                    Flag4ub=seeds(j,:)>dmax;
                    Flag4lb=seeds(j,:)<dmin;
                    seeds(j,:)=(seeds(j,:).*(~(Flag4ub+Flag4lb)))+dmax.*Flag4ub+dmin.*Flag4lb;
                end
            else
                for j=1:Number_of_seeds
                    % JJH： 新方法，强调r的提取，应该是这颗树最近的5颗树之一。
                    tnk=fix(rand*(N-nk-1))+1;%远的树（1~24）
                    while 1%保证选的是两棵不一样的树
                        rnk1=fix(rand*nk)+1;%近的树1（1~5）
                        rnk2=fix(rand*nk)+1;%近的树2
                        if (rnk1 ~=rnk2 )
                            break;
                        end
                    end
                    r1=BNN_tree.index(rnk1);%选树(1~5)
                    r2=BNN_tree.index(rnk2);
                    tn=tBNN_tree.index(tnk);%选树(1~24)
                    for d=1:D
                        if(rand<ST(i))
                            seeds(j,d)=trees(i,d)+(trees(r1,d)-trees(r2,d))*sin(rand*pi*2);%Eq.11
                            %以两棵较近树播撒种子
                        else
                            seeds(j,d)=(trees(i,d)+best_params(d))/2+(trees(tn,d)-best_params(d))*cos(rand*pi*2);%Eq.12
                            %以最佳树为指引 向比较广阔的地方播撒种子
                        end
                    end;
                    Flag4ub=seeds(j,:)>dmax;
                    Flag4lb=seeds(j,:)<dmin;
                    seeds(j,:)=(seeds(j,:).*(~(Flag4ub+Flag4lb)))+dmax.*Flag4ub+dmin.*Flag4lb;
                end;
            end;
            %%
            %Find the best one of seeds
            obj_seeds=feval(fhd,seeds',varargin{:});
            [min2num,min2num_indis]=min(obj_seeds);
            if(min2num<obj_trees(i))
                trees(i,:)=seeds(min2num_indis,:);
                obj_trees(i)=min2num;
                flin(i) = 0;%发生种子替换的奖励 我们不应该再已flin == 10 去替换这个树而是探索这个树的周边
                %种子替换之后我们应该让树在周边生成种子
                ctx(i) = 1;%3点改
            end;
            fes=fes+Number_of_seeds;
        end;
        % Update the destination if there is a better solution
        [min_tree,min_tree_index]=min(obj_trees);
        if(min_tree<best)
            best=min_tree;
            best_params=trees(min_tree_index,:);
        end;
        convergence(iter)=best;
    end;
    gbest =best_params;
    gbestval =best;
    fitcount=fes;% Read elapsed time from stopwatch
    % disp(['近邻区域比 = ',num2str(jinlin/(jinlin+fjinlin))])
    % disp(['树的迁移比 = ',num2str(qianyi/(qianyi+fqianyi))])
    time=toc;
end