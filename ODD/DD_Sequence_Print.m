function [] = DD_Sequence_Print(N, dim, n, NPerTime, exitTolerance)
%Programm that provides the best DD_Sequence by means of two different
%algorithms, DD_Sequence_max_eig and DD_Sequence_unit and writes the results in
%corresponding txt files. For each of the n times, NPertime Hamiltonians on
%a dim^2xdim^2 system are sampled and a Comb is produced for each of them,
%for which the operations are optimized.

%Sample NPerTime Hamiltonians for which the Combs are created
Hamiltonians = cell(1,NPerTime);
for i = 1:NPerTime
    Ham = rand(dim^2,dim^2) + 1j*rand(dim^2,dim^2);
    Ham = Ham + Ham';
    Ham = Ham/norm(Ham);
    Hamiltonians{i} = Ham;
end


%Time differences between slots for which to compute values
Times = logspace(-3,2,n);

%Initialize Dictionaries to store the data properly
DD_Data_eig_opt_dict = containers.Map;
DD_Data_dist_unit_opt_dict = containers.Map;
DD_Hamiltonians_dict = containers.Map;
DD_Operations_Optimized_eig_dict = containers.Map;
DD_Operations_Optimized_unitary_dict = containers.Map;
DD_Init_States_dict = containers.Map;


%Optimize and print to file
for i = 1:n
    disp(strcat('Progress: ', sprintf('%0.1f',100*(i/n)),'%'))
    DelT = Times(i);
    
    %execute operations in parallel
    for j = 1:NPerTime
        ind = int2str(NPerTime*(i-1) +j)
        Ham = Hamiltonians{i};
        [Comb, State] = Rand_comb(N, Ham, DelT);
        
        
        %Run first algorithm and write results to dictionary
        [PurityOpt, PurityStand, MutInfOpt, MutInfStand, Purity, OpEig1, OpEig2, OpEig3] = DD_Sequence_max_eig(Comb, N, dim, exitTolerance);
        Values1 = [round(str2num(ind)) DelT PurityOpt PurityStand MutInfOpt MutInfStand Purity];
        DD_Data_eig_opt_dict(ind) = Values1;
        
        %Run second algorithm and write results to dictionary
        [PurityOpt, PurityStand, MutInfOpt, MutInfStand, Purity, OpUnit1, OpUnit2, OpUnit3] = DD_Sequence_unitary(Comb, N, dim, exitTolerance);
        DD_Data_dist_unit_opt_dict(ind) = [round(str2num(ind)) DelT PurityOpt PurityStand MutInfOpt MutInfStand Purity];
        
        
        %store Hamiltonian in dictionary
        DD_Hamiltonians_dict(ind) = Ham;
        
        %store States in dictionary
        DD_Init_States_dict(ind) = State;
        
        %store operations from eigenvalue optimization in dictionary
        OpEig = cell(1,3);
        OpEig{1} = OpEig1;
        OpEig{2} = OpEig2;
        OpEig{3} = OpEig3;
        DD_Operations_Optimized_eig_dict(ind) = OpEig;
        
        
        %store operations from unitary distance optimization in dictionary
        OpUnit = cell(1,3);
        OpUnit{1} = OpUnit1;
        OpUnit{2} = OpUnit2;
        OpUnit{3} = OpUnit3;
        DD_Operations_Optimized_unitary_dict(ind) = OpUnit;
        
    end
end
     
%Set up text file for the data
Name1 = strcat('DD_Data_eig_opt','.txt');
Name2 = strcat('DD_Data_dist_unit_opt','.txt');
Name3 = strcat('DD_Hamiltonians','.txt');
Name4 = strcat('DD_Operations_Optimized_eig','.txt');
Name5 = strcat('DD_Operations_Optimized_unitary','.txt');
Name6 = strcat('DD_Init_States','.txt');

%Write headers for the files
fileID1 = fopen(Name1,'w');
fprintf(fileID1,'#Optimization obtained via Optimization of the largest Eigenvalue\n#\n', 'w');
fprintf(fileID1,'#Identifier DeltaT PurityOpt PurityStand MutInfOpt MutInfStand Purity1 Purity2 Purity3\n#\n', 'w');
fclose(fileID1);

fileID2 = fopen(Name2,'w');
fprintf(fileID2,'#Optimization obtained via Optimization of the distance to the set of unitaries\n#\n', 'w');
fprintf(fileID2,'#Identifier DeltaT PurityOpt PurityStand MutInfOpt MutInfStand Purity1 Purity2 Purity3\n#\n', 'w');
fclose(fileID2);

fileID3 = fopen(Name3,'w');
fprintf(fileID3,'#Hamiltonians used for the computation. Identifiers and Time between Slots given\n#\n', 'w');
fclose(fileID3);

fileID4 = fopen(Name4,'w');
fprintf(fileID4,'#Decoupling operations from eigenvalue optimization\n#\n', 'w');
fclose(fileID4);

fileID5 = fopen(Name5,'w');
fprintf(fileID5,'#Decoupling operations from optimizing distance to unitaries\n#\n', 'w');
fclose(fileID5);

fileID6 = fopen(Name6,'w');
fprintf(fileID6,'#Environment states used for the computation. Identifiers and Time between Slots given\n#\n', 'w');
fclose(fileID6);
        
%Iterate over the keys of the dictionaries and write content in the
%different files

for ind = 1:n*NPerTime %Array with all keys
   label = int2str(ind)
   
   %Store Values obtained from first algorithm
   Values1 = DD_Data_eig_opt_dict(label);
   DelT = Values1(2);
   fileID1 = fopen(Name1,'a');
   fprintf(fileID1, '%d %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n', Values1);
   fclose(fileID1);
   
   %Store Values obtained form second algorithm
   Values2 = DD_Data_dist_unit_opt_dict(label);
   fileID2 = fopen(Name2,'a');
   fprintf(fileID2, '%d %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n', Values2);
   fclose(fileID2);
   
   lab = strcat('[', label ,']');
   
   %Store Hamiltonians
   fileID3 = fopen(Name3,'a');        
   fprintf(fileID3, [lab ' ' num2str(DelT) '\n'], 'a');
   fclose(fileID3);
   Hamil = DD_Hamiltonians_dict(label);
   dlmwrite(Name3,Hamil,'Delimiter', ' ', '-append');
   fileID3 = fopen(Name3,'a'); 
   fprintf(fileID3, '\n', 'a');
   fclose(fileID3);
   
   %store States
   fileID6 = fopen(Name6,'a');        
   fprintf(fileID6, [lab ' ' num2str(DelT) '\n'], 'a');
   fclose(fileID6);
   State = DD_Init_States_dict(label);
   dlmwrite(Name6,State,'Delimiter', ' ', '-append');
   fileID6 = fopen(Name6,'a'); 
   fprintf(fileID6, '\n', 'a');
   fclose(fileID6);
   
   %store operations from eigenvalue optimization
   fileID4 = fopen(Name4,'a');        
   fprintf(fileID4, lab, 'a');
   fprintf(fileID4, '\n', 'a');
   fclose(fileID4);
   OpEig = DD_Operations_Optimized_eig_dict(label);
   OpEig1 = OpEig{1};
   OpEig2 = OpEig{2};
   OpEig3 = OpEig{3};
   dlmwrite(Name4,OpEig1, 'Delimiter', ' ', '-append');
   fileID4 = fopen(Name4,'a');  
   fprintf(fileID4, '\n', 'a');
   fclose(fileID4);
   dlmwrite(Name4, OpEig2, 'Delimiter', ' ','-append');
   fileID4 = fopen(Name4,'a');  
   fprintf(fileID4, '\n', 'a');
   fclose(fileID4);
   dlmwrite(Name4,OpEig3,'Delimiter', ' ', '-append');
   fileID4 = fopen(Name4,'a');
   fprintf(fileID4, '\n', 'a');
   fclose(fileID4);
        
   %store operations from unitary distance optimization
   fileID5 = fopen(Name5,'a');        
   fprintf(fileID5, lab, 'a');
   fprintf(fileID4, '\n', 'a');
   fclose(fileID5);
   OpUnit = DD_Operations_Optimized_unitary_dict(label);
   OpUnit1 = OpUnit{1};
   OpUnit2 = OpUnit{2};
   OpUnit3 = OpUnit{3};
   dlmwrite(Name5,OpUnit1,  'Delimiter', ' ','-append');
   fileID5 = fopen(Name5,'a');
   fprintf(fileID5, '\n', 'a');
   fclose(fileID5);
   dlmwrite(Name5,OpUnit2,'Delimiter', ' ', '-append');
   fileID5 = fopen(Name5,'a');
   fprintf(fileID5, '\n\n', 'a');
   fclose(fileID5);
   dlmwrite(Name5, OpUnit3, 'Delimiter', ' ', '-append');
   fileID5 = fopen(Name5,'a');
   fprintf(fileID5, '\n', 'a');
   fclose(fileID5);
end

        
        
        
        
