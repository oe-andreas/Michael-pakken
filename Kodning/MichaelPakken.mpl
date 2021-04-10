Michael := module()
options package;
export vop, prik, prikc, kryds, længde, grad, div, rot, Hessematrix, det,rum, GetJacobi, TrappeMetode, MyConstants;
prikc := (x, y) -> LinearAlgebra[DotProduct](convert(x, Vector), convert(y, Vector), conjugate = true);
prik := (x, y) -> LinearAlgebra[DotProduct](convert(x, Vector), convert(y, Vector), conjugate = false);
kryds := (x, y) -> convert(VectorCalculus[CrossProduct](convert(x, Vector), convert(y, Vector)), Vector);
rum := (x,y,z) -> prik(kryds(x,y),z);
længde := x -> LinearAlgebra[Norm](convert(x, Vector), 2);
vop := proc(X) op(convert(X, list)); end proc;
det := A -> LinearAlgebra[Determinant](A);

grad:=(X,Y)->convert(linalg[grad](X,Y),Vector[column]):
div:=V->VectorCalculus[Divergence](convert(V,Vector)):
rot:=proc(X) uses VectorCalculus;BasisFormat(false);Curl(convert(X,Vector))  end proc:

Hessematrix := proc(expr, variables := []) local vars;
vars := variables; if vars = [] then vars := indets(expr, 'name'); end if;
return VectorCalculus[Hessian](expr, convert(vars, list));
end proc;

GetJacobi := proc(parameterfremstilling, vars) local r, ru, rv, rw;
r := convert(parameterfremstilling, Vector);
if nops(vars) = 1 then
  ru := diff(r, vars[1]);
  return længde(ru);
elif nops(vars) = 2 then
  ru := diff(r, vars[1]);
  rv := diff(r, vars[2]);
  if numelems(r) = 2 then
    return abs(det(<ru | rv>));
  elif numelems(r) = 3 then
    return længde(kryds(ru, rv));
  end if;
elif nops(vars) = 3 then
  ru := diff(r, vars[1]);
  rv := diff(r, vars[2]);
  rw := diff(r, vars[3]);
  return abs(rum(ru, rv, rw));
end if;
end proc;

TrappeMetode := proc(V) local x, y, z, t;
if op(eval(V))[3] = operator then return simplify(int(V(t, 0)[1], t = 0 .. x) + int(V(x, t)[2], t = 0 .. y));
elif op(eval(V))[4] = operator then return simplify(int(V(t, 0, 0)[1], t = 0 .. x) + int(V(x, t, 0)[2], t = 0 .. y) + int(V(x, y, t)[3], t = 0 .. z));
else print("V skal være en funktion, ikke et udtryk (f.eks. skal der ikke stå V(x,y,z), men bare V)");
end if;
end proc;


MyConstants := proc(constant); #taken from appendix G of University Physics
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
else print("Du skal angive det almindelige symbol for konstanten med citationstegn")
end if;
end proc;


end module;
