<<"!math-grid info"
SetDirectory["/storage/plzen1/home/yawa/DEMES"]
MakeGamete[r_,n_Integer]:=Module[{x},
x={Random[Integer]};
Do[AppendTo[x,If[Random[]<r,1-x[[-1]],x[[-1]]]],{n-1}];
x];
MakeGamete[{h1_List,h2_List},r_]:=
h1+(h2-h1)MakeGamete[r,Length[h1]];
outer=Outer[Times,{x1,x2,x3,x4,x5,x6},{x1,x2,x3,x4,x5,x6}];
Flatten[{Diagonal[outer],2Diagonal[outer,1],2Diagonal[outer,2],2Diagonal[outer,3],2Diagonal[outer,4],2Diagonal[outer,5]}]/.{x1->0,x4->1,x3->0,x2->0,x5->0,x6->0}
gvector=RandomInteger[MultinomialDistribution[10000,%]]
uop=Join[Array[{{0,0},{0,0}}&,gvector[[1]]],Array[{{0},{0}}&,gvector[[2]]],Array[{{1,0},{1,0}}&,gvector[[3]]],Array[{{1},{1}}&,gvector[[4]]],Array[{{0,1},{0,1}}&,gvector[[5]]],Array[{{1,1},{1,1}}&,gvector[[6]]],Array[{{0,0},{0}}&,gvector[[7]]],Array[{{0},{1,0}}&,gvector[[8]]],Array[{{1,0},{1}}&,gvector[[9]]],Array[{{1},{0,1}}&,gvector[[10]]],Array[{{0,1},{1,1}}&,gvector[[11]]],Array[{{0,0},{1,0}}&,gvector[[12]]],Array[{{0},{1}}&,gvector[[13]]],Array[{{1,0},{0,1}}&,gvector[[14]]],Array[{{1},{1,1}}&,gvector[[15]]],Array[{{0,0},{1}}&,gvector[[16]]],Array[{{0},{0,1}}&,gvector[[17]]],Array[{{1,0},{1,1}}&,gvector[[18]]],Array[{{0,0},{0,1}}&,gvector[[19]]],Array[{{0},{1,1}}&,gvector[[20]]],Array[{{0,0},{1,1}}&,gvector[[21]]]];
outer=Outer[Times,{x1,x2,x3,x4,x5,x6},{x1,x2,x3,x4,x5,x6}];
Flatten[{Diagonal[outer],2Diagonal[outer,1],2Diagonal[outer,2],2Diagonal[outer,3],2Diagonal[outer,4],2Diagonal[outer,5]}]/.{x1->0,x4->0,x3->0,x2->1,x5->0,x6->0}
gvector=RandomInteger[MultinomialDistribution[10000,%]]
pop=Join[Array[{{0,0},{0,0}}&,gvector[[1]]],Array[{{0},{0}}&,gvector[[2]]],Array[{{1,0},{1,0}}&,gvector[[3]]],Array[{{1},{1}}&,gvector[[4]]],Array[{{0,1},{0,1}}&,gvector[[5]]],Array[{{1,1},{1,1}}&,gvector[[6]]],Array[{{0,0},{0}}&,gvector[[7]]],Array[{{0},{1,0}}&,gvector[[8]]],Array[{{1,0},{1}}&,gvector[[9]]],Array[{{1},{0,1}}&,gvector[[10]]],Array[{{0,1},{1,1}}&,gvector[[11]]],Array[{{0,0},{1,0}}&,gvector[[12]]],Array[{{0},{1}}&,gvector[[13]]],Array[{{1,0},{0,1}}&,gvector[[14]]],Array[{{1},{1,1}}&,gvector[[15]]],Array[{{0,0},{1}}&,gvector[[16]]],Array[{{0},{0,1}}&,gvector[[17]]],Array[{{1,0},{1,1}}&,gvector[[18]]],Array[{{0,0},{0,1}}&,gvector[[19]]],Array[{{0},{1,1}}&,gvector[[20]]],Array[{{0,0},{1,1}}&,gvector[[21]]]];
IterateDemes5[{pop_,uop_},{N1_,N2_},{m1_,m2_},{s1_,s2_},{h1_,h2_},{k1_,k2_},r_,δ_]:=Module[{p2,p3,p4,p5,u2,u3,u4,u5,rn,p6,u6},
(* choose emmigrant positions in each deme*)
rn={RandomSample[Array[{#+1}&,N1,0],Round[N1 m1]],RandomSample[Array[{#+1}&,N2,0],Round[N2 m2]]};
(* migration *)
p2=Join[Delete[pop,rn[[1]]],Part[uop,Flatten[rn[[2]]]]];
(* set aside the "pure" local genotypes (assume these never die) and the rest *)
p3=Cases[p2,{{0},{0}}|{{0,0},{0,0}}|{{0,0},{0}}|{{0},{0,0}}|{{0},{0,0,0}}|{{0,0,0},{0}}|{{0,0},{0,0,0}}|{{0,0,0},{0,0}}|{{0,0,0},{0,0,0}}];p4=Cases[p2,Except[{{0},{0}}|{{0,0},{0,0}}|{{0,0},{0}}|{{0},{0,0}}|{{0},{0,0,0}}|{{0,0,0},{0}}|{{0,0},{0,0,0}}|{{0,0,0},{0,0}}|{{0,0,0},{0,0,0}}]];

u2=Join[Delete[uop,rn[[2]]],Part[pop,Flatten[rn[[1]]]]];
u3=Cases[u2,{{1},{1}}|{{1,1},{1,1}}|{{1,1},{1}}|{{1},{1,1}}|{{1},{1,1,1}}|{{1,1,1},{1}}|{{1,1,1},{1,1}}|{{1,1},{1,1,1}}|{{1,1,1},{1,1,1}}];
u4=Cases[u2,Except[{{1},{1}}|{{1,1},{1,1}}|{{1,1},{1}}|{{1},{1,1}}|{{1},{1,1,1}}|{{1,1,1},{1}}|{{1,1,1},{1,1}}|{{1,1},{1,1,1}}|{{1,1,1},{1,1,1}}]];

Module[{l,nl,a1,a2,b1,A1,A2},(* create a True/False survival mask for all heterozygotes - perhaps this needs changing to increase speed *)
l=Table[RandomInteger[BernoulliDistribution[1-If[Count[Flatten[p4[[i]]],0]!=0,h1 s1 Count[Flatten[p4[[i]]],1]^δ Count[Flatten[p4[[i]]],0]^-k1,Count[Flatten[p4[[i]]],1]^δ s1]]]==0,{i,Length[p4]}];(* cull the unfit ind: fitness=1-s*h*(n_foreign^σ/n_local^k) or 1- s*n_foreign^σ *)
nl={};MapIndexed[If[#1,AppendTo[nl,#2]]&,l];p5=Join[p3,Delete[p4,nl]];

a1=Array[RandomInteger[{1,Length[p5]}]&,Length[p5]];A1=p5[[a1]];a2=Array[RandomInteger[{1,Length[p5]}]&,Length[p5]];A2=p5[[a2]];

p6=Transpose[{DeleteCases[MakeGamete[#,r],Null]&/@Table[If[Length[#]<3,If[Length[#]<2,Join[#,{Null,Null}],Append[#,Null]],If[Length[#]>3,DeleteCases[#,0|1,1,1],#]]&/@A1[[i]],{i,Length[A1]}],

DeleteCases[MakeGamete[#,r],Null]&/@Table[If[Length[#]<3,If[Length[#]<2,Join[#,{Null,Null}],Append[#,Null]],If[Length[#]>3,DeleteCases[#,0|1,1,1],#]]&/@A2[[i]],{i,Length[A2]}]}];
];


Module[{l,nl,a1,a2,A1,A2},
l=Table[RandomInteger[BernoulliDistribution[1-If[Count[Flatten[u4[[i]]],1]!=0,h2 s2 Count[Flatten[u4[[i]]],0]^δ Count[Flatten[u4[[i]]],1]^-k2,Count[Flatten[u4[[i]]],0]^δ s2]]]==0,{i,Length[u4]}];
nl={};MapIndexed[If[#1,AppendTo[nl,#2]]&,l];u5=Join[u3,Delete[u4,nl]];

a1=Array[RandomInteger[{1,Length[u5]}]&,Length[u5]];A1=u5[[a1]];a2=Array[RandomInteger[{1,Length[u5]}]&,Length[u5]];A2=u5[[a2]];

u6=Transpose[{DeleteCases[MakeGamete[#,r],Null]&/@Table[If[Length[#]<3,If[Length[#]<2,Join[#,{Null,Null}],Append[#,Null]],If[Length[#]>3,DeleteCases[#,1|0,1,1],#]]&/@A1[[i]],{i,Length[A1]}],

DeleteCases[MakeGamete[#,r],Null]&/@Table[If[Length[#]<3,If[Length[#]<2,Join[#,{Null,Null}],Append[#,Null]],If[Length[#]>3,DeleteCases[#,1|0,1,1],#]]&/@A2[[i]],{i,Length[A2]}]}];

];
If[Length[p6]!=N1,If[Length[p6]<N1,While[Length[p6]<N1,AppendTo[p6,p6[[RandomInteger[{1,Length[p6]}]]]]],p6=Delete[p6,Array[{RandomInteger[{1,Length[p6]-N1}]}&,Length[p6]-N1]]]];

If[Length[u6]!=N2,If[Length[u6]<N2,While[Length[u6]<N2,AppendTo[u6,u6[[RandomInteger[{1,Length[u6]}]]]]],u6=Delete[u6,Array[{RandomInteger[{1,Length[u6]-N2}]}&,Length[u6]-N2]]]];

{p6,u6}
]

IterateDemes4[{pop_, uop_}, {N1_, N2_}, {m1_, m2_}, {s1_, s2_}, {h1_, 
   h2_}, {k1_, k2_}, r_] := 
 Module[{p2, p3, p4, p5, u2, u3, u4, u5, rn, p6, u6},
  
  rn = {RandomSample[Array[{# + 1} &, N1, 0], Round[N1 m1]], 
    RandomSample[Array[{# + 1} &, N2, 0], Round[N2 m2]]};
  
  p2 = Join[Delete[pop, rn[[1]]], Part[uop, Flatten[rn[[2]]]]]; 
  p3 = Cases[
    p2, {{0}, {0}} | {{0, 0}, {0, 0}} | {{0, 0}, {0}} | {{0}, {0, 
       0}} | {{0}, {0, 0, 0}} | {{0, 0, 0}, {0}} | {{0, 0}, {0, 0, 
       0}} | {{0, 0, 0}, {0, 0}} | {{0, 0, 0}, {0, 0, 0}}]; 
  p4 = Cases[p2, 
    Except[{{0}, {0}} | {{0, 0}, {0, 0}} | {{0, 0}, {0}} | {{0}, {0, 
        0}} | {{0}, {0, 0, 0}} | {{0, 0, 0}, {0}} | {{0, 0}, {0, 0, 
        0}} | {{0, 0, 0}, {0, 0}} | {{0, 0, 0}, {0, 0, 0}}]];
  
  u2 = Join[Delete[uop, rn[[2]]], Part[pop, Flatten[rn[[1]]]]];
  u3 = Cases[
    u2, {{1}, {1}} | {{1, 1}, {1, 1}} | {{1, 1}, {1}} | {{1}, {1, 
       1}} | {{1}, {1, 1, 1}} | {{1, 1, 1}, {1}} | {{1, 1, 1}, {1, 
       1}} | {{1, 1}, {1, 1, 1}} | {{1, 1, 1}, {1, 1, 1}}];
  u4 = Cases[u2, 
    Except[{{1}, {1}} | {{1, 1}, {1, 1}} | {{1, 1}, {1}} | {{1}, {1, 
        1}} | {{1}, {1, 1, 1}} | {{1, 1, 1}, {1}} | {{1, 1, 1}, {1, 
        1}} | {{1, 1}, {1, 1, 1}} | {{1, 1, 1}, {1, 1, 1}}]];
  
  Module[{l, nl, a1, a2, b1, A1, A2},
   l = Table[
     RandomInteger[
       BernoulliDistribution[
        1 - If[Count[Flatten[p4[[i]]], 0] != 0, 
          h1 s1 Count[Flatten[p4[[i]]], 0]^-k1, s1]]] == 0, {i, 
      Length[p4]}];
   nl = {}; MapIndexed[If[#1, AppendTo[nl, #2]] &, l]; 
   p5 = Join[p3, Delete[p4, nl]];
   
   a1 = Array[RandomInteger[{1, Length[p5]}] &, Length[p5]]; 
   A1 = p5[[a1]]; 
   a2 = Array[RandomInteger[{1, Length[p5]}] &, Length[p5]]; 
   A2 = p5[[a2]];
   
   p6 = Transpose[{DeleteCases[MakeGamete[#, r], Null] & /@ 
       Table[If[Length[#] < 3, 
           If[Length[#] < 2, Join[#, {Null, Null}], 
            Append[#, Null]], #] & /@ A1[[i]], {i, Length[A1]}],
      
      DeleteCases[MakeGamete[#, r], Null] & /@ 
       Table[
        If[Length[#] < 3, 
           If[Length[#] < 2, Join[#, {Null, Null}], 
            Append[#, Null]], #] & /@ A2[[i]], {i, Length[A2]}]}];
   ];
  
  
  Module[{l, nl, a1, a2, A1, A2},
   l = Table[
     RandomInteger[
       BernoulliDistribution[
        1 - If[Count[Flatten[u4[[i]]], 1] != 0, 
          h2 s2 Count[Flatten[u4[[i]]], 1]^-k2, s2]]] == 0, {i, 
      Length[u4]}];
   nl = {}; MapIndexed[If[#1, AppendTo[nl, #2]] &, l]; 
   u5 = Join[u3, Delete[u4, nl]];
   
   a1 = Array[RandomInteger[{1, Length[u5]}] &, Length[u5]]; 
   A1 = u5[[a1]]; 
   a2 = Array[RandomInteger[{1, Length[u5]}] &, Length[u5]]; 
   A2 = u5[[a2]];
   
   u6 = Transpose[{DeleteCases[MakeGamete[#, r], Null] & /@ 
       Table[If[Length[#] < 3, 
           If[Length[#] < 2, Join[#, {Null, Null}], 
            Append[#, Null]], #] & /@ A1[[i]], {i, Length[A1]}],
      
      DeleteCases[MakeGamete[#, r], Null] & /@ 
       Table[If[Length[#] < 3, 
           If[Length[#] < 2, Join[#, {Null, Null}], 
            Append[#, Null]], #] & /@ A2[[i]], {i, Length[A2]}]}];
   
   ];
  If[Length[p6] != N1, 
   If[Length[p6] < N1, 
    While[Length[p6] < N1, 
     AppendTo[p6, p6[[RandomInteger[{1, Length[p6]}]]]]], 
    p6 = Delete[p6, 
      Array[{RandomInteger[{1, Length[p6] - N1}]} &, 
       Length[p6] - N1]]]];
  
  If[Length[u6] != N2, 
   If[Length[u6] < N2, 
    While[Length[u6] < N2, 
     AppendTo[u6, u6[[RandomInteger[{1, Length[u6]}]]]]], 
    u6 = Delete[u6, 
      Array[{RandomInteger[{1, Length[u6] - N2}]} &, 
       Length[u6] - N2]]]];
  
  {p6, u6}
  ]




(*ParallelDo[Clear[pp];
pp= {Nest[IterateDemes5[#,{2000,2000},{.01,.0125},{.12,.12},{1/2,1/2},{1.5,1.5},.03,1]&,{pop,uop},5]};
Do[

With[{N1=2000,N2=2000},pp=ReplacePart[pp,-1->{Join[DeleteCases[pp[[-1]][[1]],{{0,0},{_}}|{{_},{0,0}},1,1],

ReplacePart[#,Position[#,{0,0},2,1]->{0,0,0}]&[Cases[pp[[-1]][[1]],{{0,0},{_}}|{{_},{0,0}},1,1]]],

Join[DeleteCases[pp[[-1]][[2]],{{1,1},{_}}|{{_},{1,1}},1,1],

ReplacePart[#,Position[#,{1,1},2,1]->{1,1,1}]&[Cases[pp[[-1]][[1]],{{1,1},{_}}|{{_},{1,1}},1,1]]]}

]];
pp=Join[pp,{Nest[IterateDemes5[#,{2000,2000},{.01,.0125},{.12,.12},{1/2,1/2},{1,1},.03,1]&,pp[[-1]],25]}];
,{150}];

Flatten[{Transpose[Table[{Count[Flatten[pp[[i]][[#]],1],{0}]/(2Length[pp[[i]][[#]]]),
Count[Flatten[pp[[i]][[1]],#],{0,0}]/(2Length[pp[[i]][[#]]]),
Count[Flatten[pp[[i]][[1]],#],{1}]/(2Length[pp[[i]][[#]]]),Count[Flatten[pp[[i]][[1]],#],{1,1}]/(2Length[pp[[i]][[#]]]),
Count[Flatten[pp[[i]][[1]],#],{1,0}]/(2Length[pp[[i]][[#]]]),Count[Flatten[pp[[i]][[1]],#],{_,_,_}]/(2Length[pp[[i]][[#]]])},{i,Length[pp]}]]}&/@{1,2},1]>>>"/storage/plzen1/home/yawa/DEMES/TTripl2";
Print[{"round done",Plus@@Flatten[pp[[#]]]&/@{1,2}}],
{8}]*)


ParallelDo[Clear[pp];
 pp = {Nest[
    IterateDemes5[#, {2000, 2000}, {.012, .0125}, {.12, .12}, {1/2, 
       1/2}, {2, 2}, .03, .5] &, {pop, uop}, 20]};
 Do[pp = ReplacePart[
    pp, -1 -> {Insert[pp[[-1]][[1]], 0, 
       Position[pp[[-1]][[1]], 0, 3, 1]], 
      (*Insert[pp[[-1]][[2]], 1, Position[pp[[-1]][[2]], 1, 3, 1]]*)pp[[-1]][[2]]}];
  pp = Join[
    pp, {Nest[
      IterateDemes5[#, {2000, 2000}, {.012, .0125}, {.12, .12}, {1/2, 
         1/2}, {2, 2}, .03, .5] &, pp[[-1]], 25]}];, {200}];
 
 Flatten[{Transpose[
       Table[

{Count[Flatten[pp[[i]][[#]], 1], {0}]/(2 Length[pp[[i]][[#]]]), 
 Count[Flatten[pp[[i]][[#]], 1], {0, 0}]/(2 Length[pp[[i]][[#]]]), 
 Count[Flatten[pp[[i]][[#]], 1], {1}]/(2 Length[pp[[i]][[#]]]), 
 Count[Flatten[pp[[i]][[#]], 1], {1, 1}]/(2 Length[pp[[i]][[#]]]), 
 Count[Flatten[pp[[i]][[#]], 1], {1, 0}]/(2 Length[pp[[i]][[#]]]),
 Count[Flatten[pp[[i]][[#]], 1], {0, 1}]/(2 Length[pp[[i]][[#]]])},

 {i, Length[pp]}]]} & /@ {1,
      2}, 1] >>> "/storage/plzen1/home/yawa/DEMES/testNP3";
 Print[{"round done", Plus @@ Flatten[pp[[#]]] & /@ {1, 2}}], {16}]

CloseKernels[]

