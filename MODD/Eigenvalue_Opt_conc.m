function[OpEig, OpFinal, Res_Comb] = Eigenvalue_Opt_conc(Comb, Init_Env, N, m, dim, exitTolerance)

%function that returns the best sequence of decoupling operations for m
%concatenations of a given comb with N slots. Optimizes the maximum
%eigenvalue of the resulting channel via a see-saw.

%Initial environment state
Env_new = Init_Env;

%cell for decoupling operations
Conc_Seq = {};

%cell for building blocks of the resulting comb
Building_Blocks = cell(1,m);

for i = 1:m
    %Compute best operations and resulting environment state for single
    %round
    [Env_new, DD_Seq] = Eigenvalue_Opt_One_It(Comb, N, dim, Env_new, exitTolerance);

    %Add operations to cell
    Seq = 1; 
    for j = 1:length(DD_Seq)
        Seq = TnProduct(Seq,DD_Seq{j});
        Conc_Seq{end+1} = DD_Seq{j};
    end

    %Add resulting building block to cell
    Building_Blocks{i} = TrX(Comb*TnProduct(eye(dim^2), transpose(Seq), eye(dim^2)),[2],[dim^2 dim^(2*N) dim^2]);
end

%Contract the building blocks to build new comb
Res_Comb = Building_Blocks{1};
for k = 1:m-1
    %Dimension to pad out building block with 
    Dim_new = size(Res_Comb,1)/dim;
    
    %pad out old comb
    Res_Comb = TnProduct(Res_Comb, eye(dim^3));
    
    %New block to be attached. Systems permuted so that lines to be
    %contracted over fit.
    Build = syspermute(Building_Blocks{k+1},[2,1,3,4], [dim, dim, dim, dim]);
    %partially transpose for link product
    Build = PartTr(Build,[1],[dim dim^3]);
    %pad out
    Build = TnProduct(eye(Dim_new),Build);
    %Contract
    Res_Comb = TrX(Res_Comb*Build, [2], [Dim_new dim dim^3]);
end



%Find best operations for the resulting comb
[Env_new, DD_Seq] = Eigenvalue_Opt_One_It(Res_Comb, N, dim, Init_Env, exitTolerance);

%Sort these operations into the list of operations
Conc_Seq_Final = {};
for l = 1:m
    for k = 1:N
        Conc_Seq_Final{end+1} = Conc_Seq{(l-1)*N + k};
    end
    if l<m
        Conc_Seq_Final{end+1} = DD_Seq{l};
    end
end

OpFinal = DD_Seq;
OpEig = Conc_Seq_Final;

%Construct resulting final comb
Env_state = TnProduct(eye(dim),transpose(Init_Env),eye(dim^(2*N +2)));
Res_Comb = TrX(Res_Comb*Env_state, [2 4], [dim dim dim^(2*N+1),dim]);
end