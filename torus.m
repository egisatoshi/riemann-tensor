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

Γ2=Table[Sum[Ig[[i,j]]Γ1[[j,k,l]],{j,M}],
         {i,M},{k,M},{l,M}] //ExpandAll//Simplify;

(* Riemann curvature tensor *)

R=Table[D[Γ2[[i,j,k]],x[[l]]]
       -D[Γ2[[i,j,l]],x[[k]]]
       +Sum[Γ2[[m,j,k]]Γ2[[i,m,l]]
           -Γ2[[m,j,l]]Γ2[[i,m,k]],
            {m,M}],
        {i,M},{j,M},{k,M},{l,M}
        ] //ExpandAll//Simplify;

cR=Table[Sum[g[[i,j]]R[[j,k,l,m]],{j,M}],
         {i,M},{k,M},{l,M},{m,M}];

(* Ricci tensor *)
Ric=Table[Sum[R[[i,j,k,i]],
          {i,M}],{j,M},{k,M}] //ExpandAll//Simplify;
  
(* Scalar curvature *)
r=Sum[Ig[[i,j]]Ric[[i,j]],
      {i,M},{j,M}] //ExpandAll//Simplify;

(* Covariant derivative of Riemann curvature tensor *)
  
DR=Table[D[cR[[i,j,k,l]],x[[m]]]
        -Sum[Γ2[[n,m,i]]cR[[n,j,k,l]],{n,M}]
        -Sum[Γ2[[n,m,j]]cR[[i,n,k,l]],{n,M}]
        -Sum[Γ2[[n,m,k]]cR[[i,j,n,l]],{n,M}]
        -Sum[Γ2[[n,m,l]]cR[[i,j,k,n]],{n,M}],
         {i,M},{j,M},{k,M},{l,M},{m,M}] //ExpandAll//Simplify;
  
(* Bianchi identities *)
  
Table[DR[[i,j,k,l,m]]+DR[[i,j,l,m,k]]+DR[[i,j,m,k,l]],
      {h,M},{i,M},{j,M},{k,M},{l,M}] //ExpandAll//Simplify;
