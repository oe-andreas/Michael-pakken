Michael := module()
options package;
export vop, prik, prikc, kryds, længde, længdec, grad, div, rot, Hessematrix, det,rum, GetJacobi, TrappeMetode, MyConstants, Diagonalize, Stamvektorfelt, vsolve, vintegrate;
prikc := (x, y) -> LinearAlgebra[DotProduct](convert(x, Vector), convert(y, Vector), conjugate = true);
prik := (x, y) -> LinearAlgebra[DotProduct](convert(x, Vector), convert(y, Vector), conjugate = false);
kryds := (x, y) -> convert(VectorCalculus[CrossProduct](convert(x, Vector), convert(y, Vector)), Vector);
rum := (x,y,z) -> prik(kryds(x,y),z);
længde := x -> LinearAlgebra[Norm](convert(x, Vector), 2, conjugate=false);
længdec := x -> LinearAlgebra[Norm](convert(x, Vector), 2);
vop := proc(X) op(convert(X, list)); end proc;
det := A -> LinearAlgebra[Determinant](A);

grad:=(X,Y)->convert(linalg[grad](X,Y),Vector[column]):
div:=V->VectorCalculus[Divergence](V):
rot:=proc(X) uses VectorCalculus;BasisFormat(false);Curl(X)  end proc:


Hessematrix := proc(expr, variables := []) local vars;
vars := variables; if vars = [] then vars := indets(expr, 'name'); end if;
return VectorCalculus[Hessian](expr, convert(vars, list));
end proc;

GetJacobi := proc(parameterfremstilling, vars) local r, assumptions, variables, i, matrix, var;
r := convert(parameterfremstilling, Vector);
assumptions := 1 = 1;
variables := [];
matrix := Matrix(numelems(r), nops(vars));
for i to nops(vars) do
   try
      var := lhs(vars[i]);
      assumptions := assumptions, var <= rhs(rhs(vars[i])), lhs(rhs(vars[i])) <= var;
   catch:
      var := vars[i];
   end try;
   matrix[1 .. numelems(r), i] := diff(r, var);
end do;
if nops(vars) = 1 then return (simplify(`længde`(matrix)) assuming assumptions);
elif nops(vars) = 2 then
   if numelems(r) = 2 then return (simplify(abs(det(matrix))) assuming assumptions);
   elif numelems(r) = 3 then return (simplify(`længde`(kryds(matrix[1 .. numelems(r), 1], matrix[1 .. numelems(r), 2]))) assuming assumptions);
   end if;
elif nops(vars) = 3 then return (simplify(abs(det(matrix))) assuming assumptions);
end if;
end proc:

TrappeMetode := proc(V::procedure) local x, y, z, t;
if op(eval(V))[3] = operator then return simplify(int(V(t, 0)[1], t = 0 .. x) + int(V(x, t)[2], t = 0 .. y));
elif op(eval(V))[4] = operator then return simplify(int(V(t, 0, 0)[1], t = 0 .. x) + int(V(x, t, 0)[2], t = 0 .. y) + int(V(x, y, t)[3], t = 0 .. z));
else print("V skal være en funktion, ikke et udtryk (f.eks. skal der ikke stå V(x,y,z), men bare V)");
end if;
end proc;

Diagonalize := proc(A, unitarily := true, positive := true) local evals, evecs, Lambda, cols, i, S, dim, col, lambda1; if unitarily and `not`(Equal((HermitianTranspose(A)) . A, A . (HermitianTranspose(A)))) then error "Not unitarily diagonalizable. Use 'unitarily=false' to try non-unitary diagonalization"; end if; evals, S := Eigenvectors(A); Lambda := DiagonalMatrix(evals); dim := op(A)[2]; if `længde`(S[1 .. dim, dim]) = 0 then error "Not diagonalizable, too few linearly independent eigenvectors"; end if; if unitarily then S := Matrix(GramSchmidt([seq(S[1 .. dim, i], i = 1 .. dim)], normalized)); end if; if positive and `not`(0 <= Determinant(S)) then col := S[1 .. dim, 1]; S[1 .. dim, 1] := S[1 .. dim, 2]; S[1 .. dim, 2] := col; lambda1 := Lambda[1, 1]; Lambda[1, 1] := Lambda[2, 2]; Lambda[2, 2] := lambda1; end if; return Lambda, S; end proc:


MyConstants := proc(constant::string); #taken from appendix G of University Physics
if constant = "R" then return MyConstants("N_A")*MyConstants("k");
elif constant = "c" then return 2.99792458*10^8*Unit(('m')/('s'));
elif constant = "e" then return 1.602176634*10^(-19)*Unit('C');
elif constant = "G" then return 6.67408*10^(-11)*Unit((('N')*('m')^2)/(('kg')^2));
elif constant = "h" then return 6.62607015*10^(-34)*Unit(('J')*('S'));
elif constant = "k" then return 1.380649*10^(-23)*Unit(('J')/('K'));
elif constant = "N_A" then return 6.02214076*10^(23)*Unit(('mol')^(-1));
elif constant = "m_e" then return 9.10938356*10^(-31)*Unit(('kg'));
elif constant = "m_p" then return 1.672621898*10^(-27)*Unit(('kg'));
elif constant = "m_n" then return 1.674927471*10^(-27)*Unit(('kg'));
elif constant = "mu_0" then return 4*Pi*10^(-7)*Unit(('Wb')/('A'*'m'));
elif constant = "epsilon_0" then return 1/(MyConstants("mu_0")*MyConstants("c")^2);
elif constant = "F" then return 9.647*10^4*Unit('C');
else print("Du skal angive det almindelige symbol for konstanten med citationstegn. R, c, e, G, h, K, N_A, m_e, m_p, m_n, mu_0, epsilon_0")
end if;
end proc;

Stamvektorfelt := proc(V::procedure, verbose := true) local W, i; W := unapply(kryds(-<x, y, z>, <seq(int(u*V(u*x, u*y, u*z)[i], u = 0 .. 1), i = 1 .. 3)>), [x, y, z]); if verbose and `not`(Equal(rot(W)(x, y, z), V(x, y, z))) then print("Der findes ikke et stamvektorfeltet, men dette er output fra formlen"); end if; return W; end proc:

vsolve := proc(vligning, other_arguments:=[]) local alleqs, i; alleqs := {seq(lhs(vligning)[i] = rhs(vligning)[i], i = 1 .. numelems(rhs(vligning)))}; if other_arguments = [] then solve(alleqs); else solve(alleqs, op(other_arguments)); end if; end proc:

vintegrate := proc(v, vars) local i; <seq(int(v[i], vars), i = 1 .. numelems(v))>; end proc:

end module;
