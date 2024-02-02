classdef FONSGAIII < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation> <constrained/none>
% Nondominated sorting genetic algorithm III

%------------------------------- Reference --------------------------------
% K. Deb and H. Jain, An evolutionary many-objective optimization algorithm
% using reference-point based non-dominated sorting approach, part I:
% Solving problems with box constraints, IEEE Transactions on Evolutionary
% Computation, 2014, 18(4): 577-601.
% Yang, X., Zou, J., Yang, S., Zheng, J., Liu, Y., 2023. A fuzzy decision variables framework 
% for large-scale multiobjective optimization. IEEE Trans. Evol. Comput. 27(3), 445-459. 
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting 
            [Rate,Acc] = Algorithm.ParameterSet(0.8,0.4);
            [Z,Problem.N] = UniformPoint(Problem.N,Problem.M);
            Population    = Problem.Initialization();
            Zmin       = min(Population.objs,[],1);
            %% Optimization
            while Algorithm.NotTerminated(Population)  
                MatingPool = MatingSelection(Population.objs,Zmin);
                OffDec     = OperatorGA(Problem,Population(MatingPool).decs);
                %% FDV
                if Problem.FE/Problem.maxFE <= Rate
                    Offspring = FDVOperator(Problem,Rate,Acc,OffDec);
                else
                    Offspring = Problem.Evaluation(OffDec);
                end
                Zmin       = min(Population.objs,[],1);
                Population = EnvironmentalSelection([Population,Offspring],Problem.N,Z,Zmin);
                
            end
            Problem.print();
            
        end
    end
end