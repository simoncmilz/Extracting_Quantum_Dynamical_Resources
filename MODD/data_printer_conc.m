function [] = data_printer_conc(values)


%Set up text file for the data
%Set up text file for the data
Name1 = strcat('DD_Data_eig_opt_Conc','.txt');
Name2 = strcat('DD_Data_dist_unit_opt_Conc','.txt');
Name3 = strcat('DD_Hamiltonians_Conc','.txt');
Name4 = strcat('DD_Operations_Optimized_eig_Conc','.txt');
Name5 = strcat('DD_Operations_Optimized_unitary_Conc','.txt');
Name6 = strcat('DD_Init_States_Conc','.txt');
Name8= strcat('DD_Lindbladians_Conc','.txt');

%Check if file is from Eigenvalue or distance to unitary optimization

if values{1} == 1
    %Convert to array
    Values1 = cell2mat(values(2:9));
    % open file and write to it
    fileID1 = fopen(Name1,'a+');
    fprintf(fileID1, '%d %.6f %.6f %.6f %.6f %.6f %.6f %.6f \n', Values1);
    fprintf(fileID1, '\n');
    fclose(fileID1);
end 

if values{1} == 2
    %Convert to array
    Values2 = cell2mat(values(2:9));
    % open file and write to it
    fileID2 = fopen(Name2,'a+');
    fprintf(fileID2, '%d %.6f %.6f %.6f %.6f %.6f %.6f %.6f \n', Values2);
    fprintf(fileID2, '\n');
    fclose(fileID2);
end 

if values{1} == 3
    identifier = strcat('[',int2str(values{2}),']','\n');
    Hamil = values{3};
    % open file and write to it
    fileID3 = fopen(Name3,'a+');
    fprintf(fileID3, identifier);
    fclose(fileID3);
    dlmwrite(Name3,Hamil,'Delimiter', ' ', '-append');
    fileID3 = fopen(Name3,'a+'); 
    fprintf(fileID3, '\n');
    fclose(fileID3);
    
end

if values{1} == 4
    identifier = strcat('[',int2str(values{2}),']','\n');
    OpEig1 = values{3};
    OpEig2 = values{4};
    OpEig3 = values{5};
    OpEig4 = values{6};
    OpEig5 = values{7};
    OpEig6 = values{8};
    OpEig7 = values{9};
    OpEig8 = values{10};
    OpEig9 = values{11};
    OpEig10 = values{12};
    OpEig11 = values{13};
    OpEig12 = values{14};
    OpEig13 = values{15};
    OpEig14 = values{16};
    OpEig15 = values{17};
    % open file and write to it
    fileID4 = fopen(Name4,'a+');
    fprintf(fileID4, identifier);
    fclose(fileID4);
    dlmwrite(Name4,OpEig1, 'Delimiter', ' ', '-append');
    fileID4 = fopen(Name4,'a+');  
    fprintf(fileID4, '\n');
    fclose(fileID4);
    dlmwrite(Name4, OpEig2, 'Delimiter', ' ','-append');
    fileID4 = fopen(Name4,'a+');  
    fprintf(fileID4, '\n');
    fclose(fileID4);
    dlmwrite(Name4,OpEig3,'Delimiter', ' ', '-append');
    fileID4 = fopen(Name4,'a+');
    fprintf(fileID4, '\n');
    fclose(fileID4); 
    dlmwrite(Name4,OpEig4,'Delimiter', ' ', '-append');
    fileID4 = fopen(Name4,'a+');
    fprintf(fileID4, '\n');
    fclose(fileID4);   
     dlmwrite(Name4,OpEig5,'Delimiter', ' ', '-append');
    fileID4 = fopen(Name4,'a+');
    fprintf(fileID4, '\n');
    fclose(fileID4);   
     dlmwrite(Name4,OpEig6,'Delimiter', ' ', '-append');
    fileID4 = fopen(Name4,'a+');
    fprintf(fileID4, '\n');
    fclose(fileID4);   
     dlmwrite(Name4,OpEig7,'Delimiter', ' ', '-append');
    fileID4 = fopen(Name4,'a+');
    fprintf(fileID4, '\n');
    fclose(fileID4);   
     dlmwrite(Name4,OpEig8,'Delimiter', ' ', '-append');
    fileID4 = fopen(Name4,'a+');
    fprintf(fileID4, '\n');
    fclose(fileID4);   
     dlmwrite(Name4,OpEig9,'Delimiter', ' ', '-append');
    fileID4 = fopen(Name4,'a+');
    fprintf(fileID4, '\n');
    fclose(fileID4);   
     dlmwrite(Name4,OpEig10,'Delimiter', ' ', '-append');
    fileID4 = fopen(Name4,'a+');
    fprintf(fileID4, '\n');
    fclose(fileID4);   
     dlmwrite(Name4,OpEig11,'Delimiter', ' ', '-append');
    fileID4 = fopen(Name4,'a+');
    fprintf(fileID4, '\n');
    fclose(fileID4);   
     dlmwrite(Name4,OpEig12,'Delimiter', ' ', '-append');
    fileID4 = fopen(Name4,'a+');
    fprintf(fileID4, '\n');
    fclose(fileID4);   
     dlmwrite(Name4,OpEig13,'Delimiter', ' ', '-append');
    fileID4 = fopen(Name4,'a+');
    fprintf(fileID4, '\n');
    fclose(fileID4);   
     dlmwrite(Name4,OpEig14,'Delimiter', ' ', '-append');
    fileID4 = fopen(Name4,'a+');
    fprintf(fileID4, '\n');
    fclose(fileID4);   
     dlmwrite(Name4,OpEig15,'Delimiter', ' ', '-append');
    fileID4 = fopen(Name4,'a+');
    fprintf(fileID4, '\n');
    fclose(fileID4);   
end

if values{1} == 5
    identifier = strcat('[',int2str(values{2}),']','\n');
    OpUnit1 = values{3};
    OpUnit2 = values{4};
    OpUnit3 = values{5};
    OpUnit4 = values{6};
    OpUnit5 = values{7};
    OpUnit6 = values{8};
    OpUnit7 = values{9};
    OpUnit8 = values{10};
    OpUnit9 = values{11};
    OpUnit10 = values{12};
    OpUnit11 = values{13};
    OpUnit12 = values{14};
    OpUnit13 = values{15};
    OpUnit14 = values{16};
    OpUnit15 = values{17};
    % open file and write to it
    fileID5 = fopen(Name5,'a+');
    fprintf(fileID5, identifier);
    fclose(fileID5);
    dlmwrite(Name5,OpUnit1, 'Delimiter', ' ', '-append');
    fileID5 = fopen(Name5,'a+');  
    fprintf(fileID5, '\n');
    fclose(fileID5);
    dlmwrite(Name5, OpUnit2, 'Delimiter', ' ','-append');
    fileID5 = fopen(Name5,'a+');  
    fprintf(fileID5, '\n');
    fclose(fileID5);
    dlmwrite(Name5,OpUnit3,'Delimiter', ' ', '-append');
    fileID5 = fopen(Name5,'a+');
    fprintf(fileID5, '\n');
    fclose(fileID5); 
    dlmwrite(Name5,OpUnit4,'Delimiter', ' ', '-append');
    fileID5 = fopen(Name5,'a+');
    fprintf(fileID5, '\n');
    fclose(fileID5); 
    dlmwrite(Name5,OpUnit5,'Delimiter', ' ', '-append');
    fileID5 = fopen(Name5,'a+');
    fprintf(fileID5, '\n');
    fclose(fileID5); 
    dlmwrite(Name5,OpUnit6,'Delimiter', ' ', '-append');
    fileID5 = fopen(Name5,'a+');
    fprintf(fileID5, '\n');
    fclose(fileID5); 
    dlmwrite(Name5,OpUnit7,'Delimiter', ' ', '-append');
    fileID5 = fopen(Name5,'a+');
    fprintf(fileID5, '\n');
    fclose(fileID5); 
    dlmwrite(Name5,OpUnit8,'Delimiter', ' ', '-append');
    fileID5 = fopen(Name5,'a+');
    fprintf(fileID5, '\n');
    fclose(fileID5); 
    dlmwrite(Name5,OpUnit9,'Delimiter', ' ', '-append');
    fileID5 = fopen(Name5,'a+');
    fprintf(fileID5, '\n');
    fclose(fileID5); 
    dlmwrite(Name5,OpUnit10,'Delimiter', ' ', '-append');
    fileID5 = fopen(Name5,'a+');
    fprintf(fileID5, '\n');
    fclose(fileID5); 
    dlmwrite(Name5,OpUnit11,'Delimiter', ' ', '-append');
    fileID5 = fopen(Name5,'a+');
    fprintf(fileID5, '\n');
    fclose(fileID5); 
    dlmwrite(Name5,OpUnit12,'Delimiter', ' ', '-append');
    fileID5 = fopen(Name5,'a+');
    fprintf(fileID5, '\n');
    fclose(fileID5); 
    dlmwrite(Name5,OpUnit13,'Delimiter', ' ', '-append');
    fileID5 = fopen(Name5,'a+');
    fprintf(fileID5, '\n');
    fclose(fileID5); 
    dlmwrite(Name5,OpUnit14,'Delimiter', ' ', '-append');
    fileID5 = fopen(Name5,'a+');
    fprintf(fileID5, '\n');
    fclose(fileID5); 
    dlmwrite(Name5,OpUnit15,'Delimiter', ' ', '-append');
    fileID5 = fopen(Name5,'a+');
    fprintf(fileID5, '\n');
    fclose(fileID5); 
end

if values{1} == 6
    identifier = strcat('[',int2str(values{2}),']','\n');
    State = values{3};
    % open file and write to it
    fileID6 = fopen(Name6,'a+');
    fprintf(fileID6, identifier);
    fclose(fileID6);
    dlmwrite(Name6,State, 'Delimiter', ' ', '-append');
    fileID6 = fopen(Name6,'a+');  
    fprintf(fileID6, '\n');
    fclose(fileID6);
   
end


if values{1} == 8
    identifier = strcat('[',int2str(values{2}),']','\n');
    Lindbladian = values{3};
    % open file and write to it
    fileID8 = fopen(Name8,'a+');
    fprintf(fileID8, identifier);
    fclose(fileID8);
    dlmwrite(Name8,Lindbladian, 'Delimiter', ' ', '-append');
    fileID8 = fopen(Name8,'a+');  
    fprintf(fileID8, '\n');
    fclose(fileID8);
   
end


end

