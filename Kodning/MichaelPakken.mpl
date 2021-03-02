Michael := module()
options package;
export vop, prik, kryds, længde, gradient, Hessematrix, det;
prik := (x, y) -> VectorCalculus[DotProduct](convert(x, Vector), convert(y, Vector));
kryds := (x, y) -> convert(VectorCalculus[CrossProduct](convert(x, Vector), convert(y, Vector)), Vector);
længde := x -> LinearAlgebra[Norm](convert(x, Vector), 2);
vop := proc(X) op(convert(X, list)); end proc;
det := A -> LinearAlgebra[Determinant](A);

gradient := proc(expr, variables := []) local grad, i, vars;
vars := variables; if vars = [] then vars := indets(expr, 'name'); end if;
grad := Vector(nops(vars));
for i to nops(vars) do grad[i] := diff(expr, vars[i]); end do;
return grad;
end proc;

Hessematrix := proc(expr, variables := []) local vars;
vars := variables; if vars = [] then vars := indets(expr, 'name'); end if;
return VectorCalculus[Hessian](expr, convert(vars, list));
end proc;

end module;
