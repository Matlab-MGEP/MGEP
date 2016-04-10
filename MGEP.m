function [Population,Cost] = MGEP()
%=======================================================================%
% A  gene  expression  programming  framework for  evolutionary design  %
% of metaheuristic algorithms by Amin Rahati and Hojjat Rakhshani       %
% Programming dates: Aug to Nov 2015                                    %
% Last revised: Jan  2016   (simplified version for demo only)          %
%=======================================================================%

%=======================================================================%
%  MGEP source codes version 1.0                                        %
%                                                                       %
%  Developed in MATLAB R2009a(7.8)                                      %
%                                                                       %
%  Authors and programmers: Amin Rahati, Hojjat Rakhshani               %
%                                                                       %
%  E-Mail: a.rahati@cs.usb.ac.ir, hojjatrakhshani@gmail.com             %
%=======================================================================%


%=======================================================================%
% This demo program only implements a minimal version of MGEP, which is %
% a modified gene expression programming framework designed to evolve   %
% existing metaheuristic algorithms. This code is specially implemented % 
% to evolve differential evolution (DE) algorithm.                      %
%=======================================================================%



% Descriptions       
% Population_size        populations size of GEP and DE algorithms 
% Population             population of the problem solutions
% Cost                   cost vector of the solutions
% Problem                functional address of optimization problem 
% Max_FEs                maximum number of function evaluations
% Tol                    tolerance value
% Best_Cost              best known cost value for the considered problem
% XVmin                  vector of lower bounds 
% XVmax                  vector of upper bounds
% D                      dimensions of the solution vector
% Head                   length of the Head domain
%
% w                      vector scales weights of criterion 
%   w(1)                 denotes weight of  moment-of-inertia criterion
%   w(2)                 denotes weight of  Hamming distance criterion
%   w(3)                 denotes weight of  cost values criterion
%
% Operators              designed Search operators: We considered
%                        the following assumptions in the
%                        Head domain of each operator.
%                        if (ith gene==1) it denotes multiply 
%                        if (ith gene==2) it denotes plus 
%                        if (ith gene==3) it denotes minus 
%                        if (ith gene==4) it denotes unary minus 
%                        if (ith gene==5) it is a dummy operator 
%                        and doesn't make any change. 
%                        
% ORF                    open reading frame; please see  GEP  paper,
%                        Ferreira, Candida, and U. Gepsoft. "What is 
%                        Gene Expression Programming." (2008). 
% OP_Fitness             fitness value of the designed search operators
% Selection_rate         selection rate in GEP
% Mutation_rate          mutation rate rate in GEP
% DFlag                  display flag
% varargin               additional arguments for Cost function. 
%                        In this example,it determines function 
%                        NO in CEC 2013
%                       

%=======================================================================%
% Note:                                                                 %
% Please run the following command in Matlab window:                    %
% mex cec13_func.cpp -DWINDOWS                                          %
%=======================================================================%


Population_size=100;
Max_FEs=100000;
Best_Cost=-1200;
Min_Cost=inf;
Tol=1.0e-8;
Problem=@cec13_func;
D=10;
XVmin = -100.*ones(1,D);
XVmax =  100.*ones(1,D);
Head=100;
w=[-0.1 -0.1 0.9];
Selection_rate=0.9;
Mutation_rate=0.05;
DFlag=1;
varargin=3;


% step 1: initialize MGEP and metaheuristic algorithm
[Operators,ORF]   = Initialize_GEP(Population_size,Head);
[Population,Cost] = Initialize_DE(Population_size,Problem,D,XVmin,XVmax,varargin);
Iterations=Population_size;
while(true)
% steps 2-3: decode chromosomes and execute metaheuristic algorithm
[Population,Cost]  = EDE(Population,Cost,Problem,D,XVmin,XVmax,Operators,ORF,varargin);
Iterations=Iterations+Population_size;
Min_Cost=min(Min_Cost,min(Cost));
% step 4: check stopping criteria
    if (Iterations>=Max_FEs) || (abs(Best_Cost-Min_Cost)<Tol) 
        return;
    end
% step 5: evaluate search operators fitness    
OP_Fitness=Fitness(Population,Operators,Cost,w,Population_size);
% steps 6-8: design new search operators by gene expression programming
[Operators ORF]=GEP(OP_Fitness,Operators,ORF,Population_size,Selection_rate,Mutation_rate);     
% optional: display results
if (DFlag==1)
fprintf('@Iteration %d best fitness value is %e\n',Iterations,Min_Cost);
end
end
end
