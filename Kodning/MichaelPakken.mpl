Michael := module()
options package;
export vop, prik, prikc, kryds, længde, grad, div, rot, Hessematrix, det,rum, GetJacobi, TrappeMetode;
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

end module;
