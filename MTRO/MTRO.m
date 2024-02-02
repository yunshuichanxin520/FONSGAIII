classdef MTRO < PROBLEM
% <multi/many> <permutation> <constrained> <large/none> 
% The multi-objective traveling salesman problem
% c --- 0 --- Correlation parameter


%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        C;  % Adjacency matrix of each map
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            % Parameter setting
         
            
            obj.data.numN=10; 
            obj.data.Cap_Ts=xlsread("Maximum transfer capacity at the node.xlsx");
            obj.data.Windows=xlsread("Time window of the node.xlsx");
            obj.data.Distance=xlsread("Distance between nodes.xlsx");
       
            obj.data.Cap_Tp=xlsread("Maximum transportation capacity between nodes.xlsx");
       
            obj.data.T=obj.data.Distance; 
            obj.encoding = 1 + zeros(1,obj.D);
            obj.lower    = zeros(1,obj.D);
            obj.upper = 1 + zeros(1,obj.D);
          
            obj.data.v=[76, 60, 20];
            


            obj.data.Tmin = 30;
            obj.data.Tmax = 80;

            % -----------------------------
            obj.data.qt = [80];
            obj.data.s = [0];
          
            obj.data.q=80;
            % -----------------------------



       
            
            for i=1:length(obj.data.Cap_Tp(:,1))
         
                no1=obj.data.Cap_Tp(i,1);
           
                no2=obj.data.Cap_Tp(i,2);
            
                for j=1:3
                    if isnan(obj.data.Cap_Tp(i,2+j))
                        obj.data.Distance(i,2+j)=nan;
                    end
                end

                

           
                
                obj.data.T(i,[3,4,5])=round(obj.data.Distance(i,3:5)./obj.data.v);
                
            end
            obj.data.CT=[0,3.09,5.23; 
                3.09,0,26.62;
                5.23,26.62,0];
            obj.data.TT=[0,1,1;        
                1,0,2;
                1,2,0];
            obj.data.ET=[0,1.56,6;     
                    1.56,0,3.12;
                    6,3.12,0];
            
          
            
     
            obj.data.E0=[0.796,0.028,0.04];
             
      
            obj.data.CW=[30,50];
            
            
            obj.data.S=1;   
            obj.data.E=10;  
            
            obj.data.alpha=0.8;        
            obj.data.beta1=0.8;        
            obj.data.beta2=0.8;        
            obj.data.beta3=0.8;         
            
            obj.data.C0=[0.3,0.2,0.1];
            
        
            obj.data.weight=[1,1];
            
         
            obj.data.maxB=100000000;  
         
            obj.data.maxE=21000;  
            %%
            %%
            
            
            obj.data.numQ=1;
         

            obj.data.dim = length(obj.data.Distance(:,1))*3;

            obj.option.lb=0;
            obj.option.ub=1;
            obj.option.dim=length(obj.data.Distance(:,1))*3;
            
            % ones，
            % ones(1,5) 
            if length(obj.option.lb)==1
          
                obj.option.lb=ones(1,obj.option.dim)*obj.option.lb;
                obj.option.ub=ones(1,obj.option.dim)*obj.option.ub;
            end
            % ----------
          
            
            obj.option.fobj=@aimFcn_1; 



        end
        function Population = Initialization(obj,N)
        %Initialization - Generate multiple initial solutions.
        %
        %   P = obj.Initialization() randomly generates the decision
        %   variables of obj.N solutions and returns the SOLUTION objects.
        %
        %   P = obj.Initialization(N) generates N solutions.
        %
        %   This function is usually called at the beginning of algorithms.
        %
        %   Example:
        %       Population = Problem.Initialization()
        
            if nargin < 2
            	N = obj.N;
            end
            PopDec = zeros(N,obj.D);
            % Type   = arrayfun(@(i)find(obj.encoding==i),1:5,'UniformOutput',false);
            % if ~isempty(Type{1})        % Real variables
            %     PopDec(:,Type{1}) = unifrnd(repmat(obj.lower(Type{1}),N,1),repmat(obj.upper(Type{1}),N,1));
            % end
            % if ~isempty(Type{2})        % Integer variables
            %     PopDec(:,Type{2}) = round(unifrnd(repmat(obj.lower(Type{2}),N,1),repmat(obj.upper(Type{2}),N,1)));
            % end
            % if ~isempty(Type{3})        % Label variables
            %     PopDec(:,Type{3}) = round(unifrnd(repmat(obj.lower(Type{3}),N,1),repmat(obj.upper(Type{3}),N,1)));
            % end
            % if ~isempty(Type{4})        % Binary variables
            %     PopDec(:,Type{4}) = logical(randi([0,1],N,length(Type{4})));
            % end
            % if ~isempty(Type{5})        % Permutation variables
            %     [~,PopDec(:,Type{5})] = sort(rand(N,length(Type{5})),2);
            % end
            for i=1:N
                PopDec(i,:)=rand(size(obj.option.lb)).*(obj.option.ub-obj.option.lb)+obj.option.lb;
            end
            % disp(PopDec)
            Population = obj.Evaluation(PopDec);
        end
        function Population = Evaluation(obj,varargin)
 
            
            PopDec = varargin{1};

            [N,D]  = size(PopDec);
            % disp(PopDec)
      

       
            PopObj=ones(N,obj.M);
            PopCon =ones(N,4); 

            for i=1:N
                % x(i,:)=rand(size(obj.option.lb)).*(obj.option.ub-obj.option.lb)+obj.option.lb;
                [~, results]=obj.aimFcn_1(PopDec(i,:), obj.data);
           
                PopObj(i, :)= [results.fit1, results.fit2, results.fit3, results.fit4];
                PopCon(i, 1) = obj.data.Tmin - results.nowT;
                PopCon(i, 2) = results.nowT - obj.data.Tmax;
                PopCon(i, 3) = results.cost_C - obj.data.maxB;
                PopCon(i, 4) = results.cost_E - obj.data.maxE;
                obj.data.result(i) = results;
            end
            % obj.data.result = result;
            % disp(obj.data.result)
            Population = SOLUTION(PopDec,PopObj,PopCon,varargin{2:end});
            % Population = SOLUTION(PopDec,PopObj,PopCon, obj.data.result);
            obj.FE     = obj.FE + length(Population);
            
        end
        
        
        % function PopObj = CalObj(obj,PopDec)
        % 
        % 
        %     [N,D]  = size(PopDec);
        %     % disp(PopDec)
     
        %     PopObj=ones(N,obj.M);
        %     PopCon =ones(N,4); 
        %     for i=1:N
        %         % x(i,:)=rand(size(obj.option.lb)).*(obj.option.ub-obj.option.lb)+obj.option.lb;
        %         [~, result]=obj.aimFcn_1(PopDec(i,:), obj.data);
     
        %         PopObj(i, :)= [result.fit1, result.fit2, result.fit3];
        %         PopCon(i, 1) = obj.data.Tmin - result.nowT;
        %         PopCon(i, 2) = result.nowT - obj.data.Tmax;
        %         PopCon(i, 3) = result.fit1 - obj.data.maxB;
        %         PopCon(i, 4) = result.fit2 - obj.data.maxE;dfd
        %     end
        %     % disp(PopObj)
        % 
        % end

        

        %% Generate a point for hypervolume calculation
        %function R = GetOptimum(obj,~)
        %    R = zeros(1,obj.M) + obj.D;
        %end
    end
end
