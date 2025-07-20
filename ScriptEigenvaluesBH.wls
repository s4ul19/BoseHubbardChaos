#!/usr/bin/env wolframscript
(* ::Package:: *)

Get["/home/saul/jungkook/MatrixComponent.wl"]


n=4;
m=4;
J1=2;
U1=3;
Sym=Normal[SymmetricSectorHamiltonian[n,m,J1,U1,"Symmetric"]];
Anti=Normal[SymmetricSectorHamiltonian[n,m,J1,U1,"Antisymmetric"]];
SymEigenval=Eigenvalues[Sym];
AntiEigenval= Eigenvalues[Anti];


Export["/home/saul/jungkook/Data/
EigenSymmetric"<>"_n_"<>ToString[NumberForm[n,{Infinity,2}]]<>
"_m_"
<>ToString[NumberForm[m,{Infinity,2}]]<>"_J_"<>ToString[NumberForm[J1,{Infinity,2}]]<>"_U_"<>ToString[NumberForm[U1,{Infinity,2}]]<>".csv",SymEigenval,"CSV"];
Export["/home/saul/jungkook/Data/EigenAntiSymmetric"<>"_n_"<>ToString[NumberForm[n,{Infinity,2}]]<>
"_m_"
<>ToString[NumberForm[m,{Infinity,2}]]<>"_J_"<>ToString[NumberForm[J1,{Infinity,2}]]<>"_U_"<>ToString[NumberForm[U1,{Infinity,2}]]<>".csv",AntiEigenval,"CSV"];


Print["Finish."];
