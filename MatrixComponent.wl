(* ::Package:: *)

BeginPackage["MatrixComponent`"]
Basis::usage="N=Numero de Atomos, M=Numero de casillas";
hamiltonianBH::usage="Parametros, en orden (Numero de atomos, numero de sitios, J,U)";
DeepHamiltonianBH::usage="Parametros, en orden (Numero de atomos, numero de sitios, J,U)"

Begin["`Private`"]

listaBase[v_,n_]:=Table[If[i==1,v,0],{i,n}]; 

Indice[lista_]:=Module[{n=Length[lista]},SelectFirst[Range[n-1],AllTrue[lista[[#+1;;n-1]],(#==0)&]&,Missing["NotFound"]]]

BaseN[lista_,n_]:=Module[{k=Indice[lista],copia=lista},
copia[[k]]-=1;
If[k+1<=Length[copia],copia[[k+1]]=n-Total[copia[[1;;k]]];];
Do[copia[[i]]=0,{i,k+2,Length[copia]}];
copia]
Repetir[f_,ini_,n_]:=NestList[f,ini,n-1]
Basis[N_,M_]:=Module[{Base1=listaBase[N,M],n1=(M+N-1)!/(N!*(M-1)!)},Repetir[BaseN[#,N]&,Base1,n1]]

Tag[lista_]:= Module[{n=Length[lista]},Sum[Sqrt[100.*i+3]*lista[[i]],{i,n}]] 


(*Calcula el arreglo de tags*)
ArrayTag[lista_]:=Module[{n=Length[lista],A1={}},
Do[AppendTo[A1,Tag[lista[[i]]]],{i,n}];A1]

Posicion[valor_,lista_]:=Module[{pos=Position[lista,valor]},If[pos==={},Return[Null],First[pos]]](*Detecta si el valor dado esta en la lista y si lo esta da la posicion en la que se halla.*)




singleMatrizAij[listB_,listTag_]:=Module[{v1=Length[listB],A=AplicarOperatorAij[listB],k1=Posicion[Tag[listB],listTag]},
Table[{{Posicion[Tag[A[[i]][[1]]],listTag],k1},A[[i]][[2]]},{i,Length[A]}]]


CrearSparseArrayAij[lista_,index_]:=Module[{pares},pares={Flatten[#[[1]]],#[[2]]}&/@lista;
SparseArray[Rule@@@pares,index]]


TotalSparseAij[listaB_,listaTag_]:=Sum[CrearSparseArrayAij[singleMatrizAij[listaB[[i]],listaTag],Length[listaB]],{i,1,Length[listaB]}]




OperatorNi[lista_,entrada_]:=Module[{k=entrada,copy=lista},
{copy,copy[[k]]*(copy[[k]]-1)}]



AplicarOperatorNi[lista_]:=Module[{resultados,copia,suma},resultados=Table[OperatorNi[lista,k],{k,Length[lista]}];
copia=resultados[[1,1]];suma=Total[resultados[[All,2]]];{copia,suma}]


singleMatrizNi[listB_,listTag_]:=Module[{v1=Length[listB],A=AplicarOperatorNi[listB],k1=Posicion[Tag[listB],listTag]},
{{k1,k1},A[[2]]}]


CrearSparseArrayNi[list_,index_]:=Module[{k=Flatten[list[[1]]],I=index,A=list[[2]]},SparseArray[{k->A},{I,I}]]


TotalSparseNi[listaB_,listaTag_]:=
Sum[CrearSparseArrayNi[singleMatrizNi[listaB[[i]],listaTag],Length[listaB]],{i,1,Length[listaB]}]


HKin[Sparse_,val_]:=Module[{J=val},
-(J/2)*(ConjugateTranspose[Sparse]+Sparse)]


HInt[Sparse_,val_] :=Module[{U=val},
(U/2)*Sparse]

hamiltonianBH[n_,m_,J_,U_]:=Module[{U1=U,J1=J},
V=Basis[n,m];
TagV=Sort[ArrayTag[Basis[n,m]]];
TKin=TotalSparseAij[V,TagV] ;
TInt=TotalSparseNi[V,TagV];
H1=HKin[TKin,J1];
H2=HInt[TInt,U1];
H1+H2 ]

OperatorAij[i_,j_,lista_]:=Module[{copy=lista,valorI,valorII},If[i==j,Return[Null]];
valorI=copy[[i]];
valorII=copy[[j]];
copy[[i]]=valorI+1;
copy[[j]]=valorII-1;
If[Min[copy]<0,Return[Null]];
{copy,Sqrt[(valorI+1.)*valorII]}]


AplicarOperatorAij[lista_]:=Module[{n=Length[lista],pares},pares=Flatten[Table[{{i,i-1},{i-1,i}},{i,2,n}],1];
DeleteCases[OperatorAij[#[[1]],#[[2]],lista]&/@pares,Null]]


DeepHamiltonianBH[n_,m_,J_,U_]:=Module[{V,assoc,len,b,neighbors,v,coef,j0,rulesKin,TKin,diag,H2},(*Genera base y asociaci\[OAcute]n vector->\[IAcute]ndice*)V=Basis[n,m];
len=Length[V];
assoc=AssociationThread[V->Range[len]];
rulesKin=Reap[Do[b=V[[i]];
neighbors=AplicarOperatorAij[b];
Do[{v,coef}=neighbors[[j]];
j0=assoc[v];
Sow[{i,j0}->coef],{j,1,Length[neighbors]}],{i,1,len}]][[2,1]];
TKin=If[rulesKin==={},SparseArray[{},{len,len}],SparseArray[rulesKin,{len,len}]];
diag=(U/2)*(Total[#*(#-1)]&/@V);


H2=SparseArray[Band[{1,1}]->diag,{len,len}];

-J*TKin+H2]



End[]
EndPackage[]



