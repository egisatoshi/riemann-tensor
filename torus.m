(* Coordinates for torus *)
  
M=2;

x={θ,φ};

X={(a*Cos[θ]+b)Cos[φ],
   (a*Cos[θ]+b)Sin[φ],
   a*Sin[θ]};

(* Local coordinates *)

e=Table[D[X[[i]],x[[j]]],{i,M},{j,M}] //ExpandAll//Simplify;

(* Metric tensor *)

g=Table[e[[i]].e[[j]],{i,M},{j,M}] //ExpandAll//Simplify;

Ig=Inverse[g] //ExpandAll//Simplify;

(* Christoffel symbol of the first kind *)
Γ1=Table[D[g[[i,j]],x[[k]]]
        +D[g[[i,k]],x[[j]]]
        -D[g[[j,k]],x[[i]]],
         {i,M},{j,M},{k,M}]/2 //ExpandAll//Simplify;

(* Christoffel symbol of the second kind *)
Γ2=Table[Sum[Ig[[i,l]]Γ1[[l,j,k]],{l,M}],
         {i,M},{j,M},{k,M}]//ExpandAll//Simplify;

(* Riemann curvature tensor *)

R=Table[D[Γ2[[h,k,i]],x[[j]]]
       -D[Γ2[[h,j,i]],x[[k]]]
       +Sum[Γ2[[h,j,l]]Γ2[[l,k,i]]
           -Γ2[[h,k,l]]Γ2[[l,j,i]],
            {l,M}],
        {h,M},{i,M},{j,M},{k,M}
        ] //ExpandAll//Simplify;

cR=Table[Sum[g[[h,l]]R4[[l,i,j,k]],{l,M}],
         {h,M},{i,M},{j,M},{k,M} ];

(* Ricci tensor *)
Ric=Table[Sum[R[[h,i,j,h]],
          {h,M}],{i,M},{j,M}] //ExpandAll//Simplify;
  
(* Scalar curvature *)
R=Sum[Ig[[i,j]]Ric[[i,j]],
      {i,M},{j,M}] //ExpandAll//Simplify;

(* Covariant derivative of Riemann curvature tensor *)
  
DR=Table[D[R[[h,i,j,k]],x[[l]]]
        +Sum[Γ2[[h,l,m]]R[[m,i,j,k]],{m,M}]
        -Sum[Γ2[[m,l,i]]R[[h,m,j,k]],{m,M}]
        -Sum[Γ2[[m,l,j]]R[[h,i,m,k]],{m,M}]
        -Sum[Γ2[[m,l,k]]R[[h,i,j,m]],{m,M}],
         {h,M},{i,M},{j,M},{k,M},{l,M} 
         ] //ExpandAll //Simplify;
  
(* Bianchi identities *)
  
Table[DR[[i,j,k,l,m]]+DR[[i,j,l,m,k]]+DR[[i,j,m,k,l]],
      {h,M},{i,M},{j,M},{k,M},{l,M}] //ExpandAll //Simplify;
