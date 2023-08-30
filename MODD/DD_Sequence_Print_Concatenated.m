function [] = DD_Sequence_Print_Concatenated(m, N, dim, n, NPerTime, exitTolerance)
%Programm that produces data for MODD and prints it to the respective txt
%files.

%Fix Solver
cvx_solver Mosek

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setting up files to be printed to
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Set up text file for the data
Name1 = strcat('DD_Data_eig_opt_Conc','.txt');
Name2 = strcat('DD_Data_dist_unit_opt_Conc','.txt');
Name3 = strcat('DD_Hamiltonians_Conc','.txt');
Name4 = strcat('DD_Operations_Optimized_eig_Conc','.txt');
Name5 = strcat('DD_Operations_Optimized_unitary_Conc','.txt');
Name6 = strcat('DD_Init_States_Conc','.txt');
Name8 = strcat('DD_Lindbladians_Conc','.txt');

fileID1 = fopen(Name1,'w');
fprintf(fileID1,'#Optimization obtained via Optimization of the largest Eigenvalue\n#\n', 'w');
fprintf(fileID1,'#Identifier DeltaT PurityOpt PurityStand PurityBase MutInfOpt MutInfStand MutInfBase\n#\n', 'w');
fclose(fileID1);

fileID2 = fopen(Name2,'w');
fprintf(fileID2,'#Optimization obtained via Optimization of the distance to the set of unitaries\n#\n', 'w');
fprintf(fileID2,'#Identifier DeltaT PurityOpt PurityStand PurityBase MutInfOpt MutInfStand MutInfBase \n#\n', 'w');
fclose(fileID2);

fileID3 = fopen(Name3,'w');
fprintf(fileID3,'#Hamiltonians used for the computation. Identifiers for the individual runs given\n#\n', 'w');
fclose(fileID3);

fileID4 = fopen(Name4,'w');
fprintf(fileID4,'#Decoupling operations from eigenvalue optimization\n#\n', 'w');
fclose(fileID4);

fileID5 = fopen(Name5,'w');
fprintf(fileID5,'#Decoupling operations from optimizing distance to unitaries\n#\n', 'w');
fclose(fileID5);

fileID6 = fopen(Name6,'w');
fprintf(fileID6,'#Environment states used for the computation. Identifiers for the run given\n#\n', 'w');
fclose(fileID6);


fileID8 = fopen(Name8,'w');
fprintf(fileID8,'#Lindbladian in vectorized form for the respective runs. Identifiers for the run given\n#\n', 'w');
fclose(fileID8);



% setup data pipeline for worker processes
pipeline = parallel.pool.DataQueue;
% define function handle that sits on end of pipeline
afterEach(pipeline, @data_printer_conc);

%Times at which to compute
Times = logspace(-2,2,n);

%Hamiltonians for which to compute the Combs
Hams = cell(1,NPerTime);
for k = 1:NPerTime
    Hamil = rand(dim^2) + 1j*rand(dim^2);
    Hamil = Hamil + Hamil';
    Hams{k} = Hamil/norm(Hamil);
end


%Optimize and print to file
for i = 1:n
    disp(strcat('Progress: ', sprintf('%0.1f',100*(i/n)),'%'))
    DelT = Times(i);
    %Initial State for the construction
    Init_Env = Rand_state(dim);
    
    parfor j = 1:NPerTime
        Hamil = Hams{j};
        ind = NPerTime*(i-1) +j;
                
        %Run algorighm to obtain Data (opt eigenvalues)
        [MutInfOptEig, PurityOptEig, OpEig, Lindblad]  = Concatenated_DD(m, N,Hamil,Init_Env, DelT, exitTolerance, 'eigenvalue');
        %Get Data for standard sequence and do-nothing sequence
        [PurityStand, MutInfStand] = StandardSequence_Conc(Lindblad, Init_Env, DelT);
        [PurityBase, MutInfBase] = BaselineSequence_Conc(Lindblad, Init_Env, DelT);
        
        %Results from optimizing largest eigenvalue
        Values1 = cell(1,9);
        Values1{1} = 1;
        Values1{2} = ind;
        Values1{3} = DelT;
        Values1{4} = PurityOptEig;
        Values1{5} = PurityStand;
        Values1{6} = PurityBase;
        Values1{7} = MutInfOptEig;
        Values1{8} = MutInfStand;
        Values1{9} = MutInfBase;
        
        %Results from optimizing distance to unitaries [PurityOpt, MutInfOpt, Purity, OpUnit1, OpUnit2, OpUnit3]
        [MutInfOptUnit, PurityOptUnit, OpUnit, Lindblad] = Concatenated_DD(m, N,Hamil,Init_Env, DelT, exitTolerance, 'dist_unit');
        Values2 = cell(1,9);
        Values2{1} = 2;
        Values2{2} = ind;
        Values2{3} = DelT;
        Values2{4} = PurityOptUnit;
        Values2{5} = PurityStand;
        Values2{6} = PurityBase;
        Values2{7} = MutInfOptUnit;
        Values2{8} = MutInfStand;
        Values2{9} = MutInfBase;
        
        %Feed forward the used Hamiltonian
        Values3 = cell(1,3);
        Values3{1} = 3;
        Values3{2} = ind;
        Values3{3} = Hamil;
    
        %Feed forward ideal operations from eigenvalue optimization
        Values4 = cell(1,17);
        Values4{1} = 4;
        Values4{2} = ind;
        for k = 1:15
            Values4{k+2} = OpEig{k};
        end
                
        %Feed forward ideal operations from distance to unitaries
        Values5 = cell(1,17);
        Values5{1} = 5;
        Values5{2} = ind;
        for k = 1:15
            Values5{k+2} = OpUnit{k};
        end
        
        %Feed forward employed Environment state
        Values6 = cell(1,3);
        Values6{1} = 6;
        Values6{2} = ind;
        Values6{3} = Init_Env;
        
        
    
        
        %Feed forward respective Lindbladian
        Values8 = cell(1,3);
        Values8{1} = 8;
        Values8{2} = ind;
        Values8{3} = Lindblad;   
 
        
        
        % Return data to multiprocessing pipeline
        send(pipeline, Values1);
        send(pipeline, Values2);
        send(pipeline, Values3);
        send(pipeline, Values4);
        send(pipeline, Values5);
        send(pipeline, Values6);
        send(pipeline, Values8);
        
    end

end
end
