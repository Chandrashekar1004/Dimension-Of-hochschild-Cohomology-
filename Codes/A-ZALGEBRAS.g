type := "AlgebraA";

# ============================================================
if type = "AlgebraA" then
    # h is variable, K = Rationals
    # x1->y, x2->z, x3->x, x4->w

    K := FunctionField(Rationals, ["h"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    h := ind[1];;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*z2 - z1*y2;
    rels[2]  := w1*x2 - x1*w2 - x1*x2;
    rels[3]  := x1*y2 - h*y1*x2;
    rels[4]  := x1*z2 - h*(y1*w2 + z1*x2);
    rels[5]  := w1*y2 - h*y1*w2;
    rels[6]  := w1*z2 + h*(2*z1*x2 + y1*w2 - z1*w2);

    rels[7]  := y2*z3 - z2*y3;
    rels[8]  := w2*x3 - x2*w3 - x2*x3;
    rels[9]  := x2*y3 - h*y2*x3;
    rels[10] := x2*z3 - h*(y2*w3 + z2*x3);
    rels[11] := w2*y3 - h*y2*w3;
    rels[12] := w2*z3 + h*(2*z2*x3 + y2*w3 - z2*w3);

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraB" then
    # p = E(4) = i fixed, h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(GaussianRationals, ["h"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    h := ind[1];;
    p := E(4) * One(K);;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 - p*x1*y2;
    rels[2]  := w1*z2 - p*z1*w2;
    rels[3]  := z1*x2 + h*(-y1*w2);
    rels[4]  := z1*y2 + h*(-x1*w2);
    rels[5]  := w1*x2 + h*(y1*z2);
    rels[6]  := w1*y2 + h*(-x1*z2);

    rels[7]  := y2*x3 - p*x2*y3;
    rels[8]  := w2*z3 - p*z2*w3;
    rels[9]  := z2*x3 + h*(-y2*w3);
    rels[10] := z2*y3 + h*(-x2*w3);
    rels[11] := w2*x3 + h*(y2*z3);
    rels[12] := w2*y3 + h*(-x2*z3);

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraC" then
    # p = E(3) fixed, h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(GaussianRationals, ["h"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    h := ind[1];;
    p := E(4) * One(K);;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := w1*z2 - p*z1*w2;
    rels[2]  := y1*x2 - p*x1*y2;
    rels[3]  := z1*x2 + h*(x1*z2 - (p^2)*y1*z2 - x1*w2 + p*y1*w2);
    rels[4]  := z1*y2 + h*(p*x1*z2 - y1*z2 - x1*w2 + p*y1*w2);
    rels[5]  := w1*x2 + h*(p*x1*z2 + 2*(p^2)*y1*z2 - p*x1*w2 + p*y1*w2);
    rels[6]  := w1*y2 + h*(p*x1*z2 - (p^2)*y1*z2 - x1*w2 + y1*w2);

    rels[7]  := w2*z3 - p*z2*w3;
    rels[8]  := y2*x3 - p*x2*y3;
    rels[9]  := z2*x3 + h*(x2*z3 - (p^2)*y2*z3 - x2*w3 + p*y2*w3);
    rels[10] := z2*y3 + h*(p*x2*z3 - y2*z3 - x2*w3 + p*y2*w3);
    rels[11] := w2*x3 + h*(p*x2*z3 + 2*(p^2)*y2*z3 - p*x2*w3 + p*y2*w3);
    rels[12] := w2*y3 + h*(p*x2*z3 - (p^2)*y2*z3 - x2*w3 + y2*w3);

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraD" then
    # p and h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(Rationals, ["p","h"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    p := ind[1];; h := ind[2];;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := x1*y2 + y1*x2;
    rels[2]  := w1*z2 - p*z1*w2;
    rels[3]  := z1*x2 + h*(p*x1*z2);
    rels[4]  := z1*y2 + h*((p*p)*y1*z2 - x1*w2);
    rels[5]  := w1*x2 - h*(p*x1*w2);
    rels[6]  := w1*y2 + h*(-x1*z2 - y1*w2);

    rels[7]  := x2*y3 + y2*x3;
    rels[8]  := w2*z3 - p*z2*w3;
    rels[9]  := z2*x3 + h*(p*x2*z3);
    rels[10] := z2*y3 + h*((p*p)*y2*z3 - x2*w3);
    rels[11] := w2*x3 - h*(p*x2*w3);
    rels[12] := w2*y3 + h*(-x2*z3 - y2*w3);

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraE" then
    # p = E(4), h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(GaussianRationals, ["h"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    h := ind[1];;
    p := E(4) * One(K);;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 + x1*y2;
    rels[2]  := w1*z2 - p*z1*w2;
    rels[3]  := z1*x2 + h*(-x1*w2 - y1*w2);
    rels[4]  := z1*y2 + h*(-x1*w2 + y1*w2);
    rels[5]  := w1*x2 + h*(x1*z2 - y1*z2);
    rels[6]  := w1*y2 + h*(-x1*z2 - y1*z2);

    rels[7]  := y2*x3 + x2*y3;
    rels[8]  := w2*z3 - p*z2*w3;
    rels[9]  := z2*x3 + h*(-x2*w3 - y2*w3);
    rels[10] := z2*y3 + h*(-x2*w3 + y2*w3);
    rels[11] := w2*x3 + h*(x2*z3 - y2*z3);
    rels[12] := w2*y3 + h*(-x2*z3 - y2*z3);

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraF" then
    # p = E(4), h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(GaussianRationals, ["h"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    h := ind[1];;
    p := E(4) * One(K);;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 + x1*y2;
    rels[2]  := w1*z2 - p*z1*w2;
    rels[3]  := z1*x2 + h*(x1*z2 + p*y1*z2 - x1*w2 + y1*w2);
    rels[4]  := z1*y2 + h*(p*x1*z2 - y1*z2 - x1*w2 - y1*w2);
    rels[5]  := w1*x2 + h*(p*x1*z2 - p*y1*z2 - p*x1*w2 - y1*w2);
    rels[6]  := w1*y2 + h*(p*x1*z2 + p*y1*z2 - x1*w2 + p*y1*w2);

    rels[7]  := y2*x3 + x2*y3;
    rels[8]  := w2*z3 - p*z2*w3;
    rels[9]  := z2*x3 + h*(x2*z3 + p*y2*z3 - x2*w3 + y2*w3);
    rels[10] := z2*y3 + h*(p*x2*z3 - y2*z3 - x2*w3 - y2*w3);
    rels[11] := w2*x3 + h*(p*x2*z3 - p*y2*z3 - p*x2*w3 - y2*w3);
    rels[12] := w2*y3 + h*(p*x2*z3 + p*y2*z3 - x2*w3 + p*y2*w3);

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraG" then
    # p, f, h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(Rationals, ["p"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    p := ind[1];;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    f:=13/27;;
    h:=1/3;;
    rels := [];

    rels[1]  := y1*x2 - x1*y2;
    rels[2]  := w1*z2 - p*z1*w2;
    rels[3]  := z1*x2 + h*(-p*x1*z2);
    rels[4]  := z1*y2 + h*(-p*x1*z2 - (p*p)*y1*z2 - x1*w2);
    rels[5]  := w1*x2 + h*(-p*x1*w2);
    rels[6]  := w1*y2 + h*(-f*x1*z2 + x1*w2 - y1*w2);

    rels[7]  := y2*x3 - x2*y3;
    rels[8]  := w2*z3 - p*z2*w3;
    rels[9]  := z2*x3 + h*(-p*x2*z3);
    rels[10] := z2*y3 + h*(-p*x2*z3 - (p*p)*y2*z3 - x2*w3);
    rels[11] := w2*x3 + h*(-p*x2*w3);
    rels[12] := w2*y3 + h*(-f*x2*z3 + x2*w3 - y2*w3);

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraH" then
    # f, h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(Rationals, ["f","h"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    f := ind[1];; h := ind[2];;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 - x1*y2 - x1*x2;
    rels[2]  := w1*z2 + z1*w2;
    rels[3]  := z1*x2 - h*x1*w2;
    rels[4]  := z1*y2 - h*f*x1*w2 - h*y1*w2;
    rels[5]  := w1*x2 - h*x1*z2;
    rels[6]  := w1*y2 - h*f*x1*z2 - h*y1*z2;

    rels[7]  := y2*x3 - x2*y3 - x2*x3;
    rels[8]  := w2*z3 + z2*w3;
    rels[9]  := z2*x3 - h*x2*w3;
    rels[10] := z2*y3 - h*f*x2*w3 - h*y2*w3;
    rels[11] := w2*x3 - h*x2*z3;
    rels[12] := w2*y3 - h*f*x2*z3 - h*y2*z3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraI" then
    # q = E(4), h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(GaussianRationals, ["h"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    h := ind[1];;
    q := E(4) * One(K);;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 - q*x1*y2;
    rels[2]  := w1*z2 + z1*w2;
    rels[3]  := z1*x2 + h*(q*x1*z2 + q*y1*z2 - x1*w2 + q*y1*w2);
    rels[4]  := z1*y2 + h*(-x1*z2 - y1*z2 - x1*w2 + q*y1*w2);
    rels[5]  := w1*x2 + h*(-x1*z2 - q*y1*z2 - q*x1*w2 + q*y1*w2);
    rels[6]  := w1*y2 + h*(x1*z2 + q*y1*z2 - x1*w2 + y1*w2);

    rels[7]  := y2*x3 - q*x2*y3;
    rels[8]  := w2*z3 + z2*w3;
    rels[9]  := z2*x3 + h*(q*x2*z3 + q*y2*z3 - x2*w3 + q*y2*w3);
    rels[10] := z2*y3 + h*(-x2*z3 - y2*z3 - x2*w3 + q*y2*w3);
    rels[11] := w2*x3 + h*(-x2*z3 - q*y2*z3 - q*x2*w3 + q*y2*w3);
    rels[12] := w2*y3 + h*(x2*z3 + q*y2*z3 - x2*w3 + y2*w3);

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraJ" then
    # q = E(4), h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(GaussianRationals, ["h"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    h := ind[1];;
    q := E(4) * One(K);;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 - q*x1*y2;
    rels[2]  := w1*z2 + z1*w2;
    rels[3]  := z1*x2 + h*(-y1*z2 - y1*w2);
    rels[4]  := z1*y2 + h*(x1*z2 - x1*w2);
    rels[5]  := w1*x2 + h*(-y1*z2 + y1*w2);
    rels[6]  := w1*y2 + h*(-x1*z2 - x1*w2);

    rels[7]  := y2*x3 - q*x2*y3;
    rels[8]  := w2*z3 + z2*w3;
    rels[9]  := z2*x3 + h*(-y2*z3 - y2*w3);
    rels[10] := z2*y3 + h*(x2*z3 - x2*w3);
    rels[11] := w2*x3 + h*(-y2*z3 + y2*w3);
    rels[12] := w2*y3 + h*(-x2*z3 - x2*w3);

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraK" then
    # q, f, h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(Rationals, ["q","f","h"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    q := ind[1];; f := ind[2];; h := ind[3];;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 - q*x1*y2;
    rels[2]  := w1*z2 + z1*w2;
    rels[3]  := z1*x2 + h*(-x1*z2);
    rels[4]  := z1*y2 + h*(-y1*w2);
    rels[5]  := w1*x2 + h*(-x1*w2);
    rels[6]  := w1*y2 + h*(-f*y1*z2);

    rels[7]  := y2*x3 - q*x2*y3;
    rels[8]  := w2*z3 + z2*w3;
    rels[9]  := z2*x3 + h*(-x2*z3);
    rels[10] := z2*y3 + h*(-y2*w3);
    rels[11] := w2*x3 + h*(-x2*w3);
    rels[12] := w2*y3 + h*(-f*y2*z3);

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraL" then
    # q, f, h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(Rationals, ["q","f","h"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    q := ind[1];; f := ind[2];; h := ind[3];;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 - q*x1*y2;
    rels[2]  := w1*z2 + z1*w2;
    rels[3]  := z1*x2 + h*(-f*x1*w2);
    rels[4]  := z1*y2 + h*(-y1*w2);
    rels[5]  := w1*x2 + h*(-f*x1*z2);
    rels[6]  := w1*y2 + h*(-y1*z2);

    rels[7]  := y2*x3 - q*x2*y3;
    rels[8]  := w2*z3 + z2*w3;
    rels[9]  := z2*x3 + h*(-f*x2*w3);
    rels[10] := z2*y3 + h*(-y2*w3);
    rels[11] := w2*x3 + h*(-f*x2*z3);
    rels[12] := w2*y3 + h*(-y2*z3);

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraM" then
    # f, h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(Rationals, ["f","h"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    f := ind[1];; h := ind[2];;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 + x1*y2;
    rels[2]  := w1*z2 + z1*w2;
    rels[3]  := z1*x2 + h*(-y1*z2 - x1*w2);
    rels[4]  := z1*y2 + h*(-f*x1*z2 + y1*w2);
    rels[5]  := w1*x2 + h*(-x1*z2 + y1*w2);
    rels[6]  := w1*y2 + h*(y1*z2 + f*x1*w2);

    rels[7]  := y2*x3 + x2*y3;
    rels[8]  := w2*z3 + z2*w3;
    rels[9]  := z2*x3 + h*(-y2*z3 - x2*w3);
    rels[10] := z2*y3 + h*(-f*x2*z3 + y2*w3);
    rels[11] := w2*x3 + h*(-x2*z3 + y2*w3);
    rels[12] := w2*y3 + h*(y2*z3 + f*x2*w3);

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraN" then
    # f, g, h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(Rationals, ["f","g","h"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    f := ind[1];; g := ind[2];; h := ind[3];;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 + x1*y2;
    rels[2]  := w1*z2 + z1*w2;
    rels[3]  := z1*x2 + h*(g*y1*z2 - f*y1*w2);
    rels[4]  := z1*y2 + h*(-g*x1*z2 - f*x1*w2);
    rels[5]  := w1*x2 + h*(-f*y1*z2 + g*y1*w2);
    rels[6]  := w1*y2 + h*(-f*x1*z2 - g*x1*w2);

    rels[7]  := y2*x3 + x2*y3;
    rels[8]  := w2*z3 + z2*w3;
    rels[9]  := z2*x3 + h*(g*y2*z3 - f*y2*w3);
    rels[10] := z2*y3 + h*(-g*x2*z3 - f*x2*w3);
    rels[11] := w2*x3 + h*(-f*y2*z3 + g*y2*w3);
    rels[12] := w2*y3 + h*(-f*x2*z3 - g*x2*w3);

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraO" then
    # f, h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(Rationals, ["f"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    f := ind[1];; 
    h :=17/31;; 

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 + x1*y2;
    rels[2]  := w1*z2 + z1*w2;
    rels[3]  := z1*x2 + h*(-x1*z2 - f*y1*w2);
    rels[4]  := z1*y2 + h*(y1*z2 - x1*w2);
    rels[5]  := w1*x2 + h*(-f*y1*z2 + x1*w2);
    rels[6]  := w1*y2 + h*(-x1*z2 - y1*w2);

    rels[7]  := y2*x3 + x2*y3;
    rels[8]  := w2*z3 + z2*w3;
    rels[9]  := z2*x3 + h*(-x2*z3 - f*y2*w3);
    rels[10] := z2*y3 + h*(y2*z3 - x2*w3);
    rels[11] := w2*x3 + h*(-f*y2*z3 + x2*w3);
    rels[12] := w2*y3 + h*(-x2*z3 - y2*w3);

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraP" then
    # f, h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(Rationals, ["f","h"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    f := ind[1];; h := ind[2];;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 + x1*y2;
    rels[2]  := w1*z2 + z1*w2;
    rels[3]  := z1*x2 + h*(-x1*w2 - f*y1*w2);
    rels[4]  := z1*y2 + h*(-x1*w2 - y1*w2);
    rels[5]  := w1*x2 + h*(-x1*z2 + f*y1*z2);
    rels[6]  := w1*y2 + h*(x1*z2 - y1*z2);

    rels[7]  := y2*x3 + x2*y3;
    rels[8]  := w2*z3 + z2*w3;
    rels[9]  := z2*x3 + h*(-x2*w3 - f*y2*w3);
    rels[10] := z2*y3 + h*(-x2*w3 - y2*w3);
    rels[11] := w2*x3 + h*(-x2*z3 + f*y2*z3);
    rels[12] := w2*y3 + h*(x2*z3 - y2*z3);

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraQ" then
    # h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(Rationals, ["h"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    h := ind[1];;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 + x1*y2;
    rels[2]  := w1*z2 + z1*w2;
    rels[3]  := z1*x2 + h*(-x1*w2);
    rels[4]  := z1*y2 + h*(-x1*z2 - y1*z2 - x1*w2);
    rels[5]  := w1*x2 + h*(x1*z2);
    rels[6]  := w1*y2 + h*(-x1*z2 + x1*w2 - y1*w2);

    rels[7]  := y2*x3 + x2*y3;
    rels[8]  := w2*z3 + z2*w3;
    rels[9]  := z2*x3 + h*(-x2*w3);
    rels[10] := z2*y3 + h*(-x2*z3 - y2*z3 - x2*w3);
    rels[11] := w2*x3 + h*(x2*z3);
    rels[12] := w2*y3 + h*(-x2*z3 + x2*w3 - y2*w3);

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraR" then
    # h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(Rationals, ["h"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    h := ind[1];;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 + x1*y2;
    rels[2]  := w1*z2 + z1*w2;
    rels[3]  := z1*x2 + h*(-x1*z2 - y1*z2 - x1*w2);
    rels[4]  := z1*y2 + h*(-x1*w2);
    rels[5  ]  := w1*x2 + h*(-y1*z2);
    rels[6]  := w1*y2 + h*(y1*z2 + x1*w2 - y1*w2);

    rels[7]  := y2*x3 + x2*y3;
    rels[8]  := w2*z3 + z2*w3;
    rels[9]  := z2*x3 + h*(-x2*z3 - y2*z3 - x2*w3);
    rels[10] := z2*y3 + h*(-x2*w3);
    rels[11] := w2*x3 + h*(-y2*z3);
    rels[12] := w2*y3 + h*(y2*z3 + x2*w3 - y2*w3);

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraS" then
    # h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(Rationals, ["h"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    h := ind[1];;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 + x1*y2;
    rels[2]  := w1*z2 + z1*w2;
    rels[3]  := z1*x2 + h*(x1*z2 - y1*z2 - x1*w2 - y1*w2);
    rels[4]  := z1*y2 + h*(-x1*z2 + y1*z2 - x1*w2 - y1*w2);
    rels[5]  := w1*x2 + h*(-x1*z2 - y1*z2 + x1*w2 - y1*w2);
    rels[6]  := w1*y2 + h*(-x1*z2 - y1*z2 - x1*w2 + y1*w2);

    rels[7]  := y2*x3 + x2*y3;
    rels[8]  := w2*z3 + z2*w3;
    rels[9]  := z2*x3 + h*(x2*z3 - y2*z3 - x2*w3 - y2*w3);
    rels[10] := z2*y3 + h*(-x2*z3 + y2*z3 - x2*w3 - y2*w3);
    rels[11] := w2*x3 + h*(-x2*z3 - y2*z3 + x2*w3 - y2*w3);
    rels[12] := w2*y3 + h*(-x2*z3 - y2*z3 - x2*w3 + y2*w3);

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraT" then
    # h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(Rationals, ["h"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    h := ind[1];;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 + x1*y2;
    rels[2]  := w1*z2 + z1*w2;
    rels[3]  := z1*x2 + h*(x1*z2 - y1*z2 - x1*w2 - y1*w2);
    rels[4]  := z1*y2 + h*(-x1*z2 + y1*z2 - x1*w2 - y1*w2);
    rels[5]  := w1*x2 + h*(-x1*z2 - y1*z2 - x1*w2 + y1*w2);
    rels[6]  := w1*y2 + h*(-x1*z2 - y1*z2 + x1*w2 - y1*w2);

    rels[7]  := y2*x3 + x2*y3;
    rels[8]  := w2*z3 + z2*w3;
    rels[9]  := z2*x3 + h*(x2*z3 - y2*z3 - x2*w3 - y2*w3);
    rels[10] := z2*y3 + h*(-x2*z3 + y2*z3 - x2*w3 - y2*w3);
    rels[11] := w2*x3 + h*(-x2*z3 - y2*z3 - x2*w3 + y2*w3);
    rels[12] := w2*y3 + h*(-x2*z3 - y2*z3 + x2*w3 - y2*w3);

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraU" then
    # h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(Rationals, ["h"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    h := ind[1];;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 + x1*y2;
    rels[2]  := w1*z2 + z1*w2;
    rels[3]  := z1*x2 + h*(x1*z2 - y1*z2 - x1*w2 - y1*w2);
    rels[4]  := z1*y2 + h*(-x1*z2 - y1*z2 - x1*w2 + y1*w2);
    rels[5]  := w1*x2 + h*(-x1*z2 - y1*z2 + x1*w2 - y1*w2);
    rels[6]  := w1*y2 + h*(-x1*z2 + y1*z2 - x1*w2 - y1*w2);

    rels[7]  := y2*x3 + x2*y3;
    rels[8]  := w2*z3 + z2*w3;
    rels[9]  := z2*x3 + h*(x2*z3 - y2*z3 - x2*w3 - y2*w3);
    rels[10] := z2*y3 + h*(-x2*z3 - y2*z3 - x2*w3 + y2*w3);
    rels[11] := w2*x3 + h*(-x2*z3 - y2*z3 + x2*w3 - y2*w3);
    rels[12] := w2*y3 + h*(-x2*z3 + y2*z3 - x2*w3 - y2*w3);

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraV" then
    # h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(Rationals, ["h"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    h := ind[1];;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 - x1*y2;
    rels[2]  := w1*z2 + z1*w2;
    rels[3]  := z1*x2 + h*(-y1*z2 - x1*w2);
    rels[4]  := z1*y2 - h*y1*z2;
    rels[5]  := w1*x2 + h*(x1*z2 - y1*z2);
    rels[6]  := w1*y2 - h*y1*w2;

    rels[7]  := y2*x3 - x2*y3;
    rels[8]  := w2*z3 + z2*w3;
    rels[9]  := z2*x3 + h*(-y2*z3 - x2*w3);
    rels[10] := z2*y3 - h*y2*z3;
    rels[11] := w2*x3 + h*(x2*z3 - y2*z3);
    rels[12] := w2*y3 - h*y2*w3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraW" then
    # f, h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(Rationals, ["f","h"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    f := ind[1];; h := ind[2];;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 - x1*y2;
    rels[2]  := w1*z2 + z1*w2;
    rels[3]  := z1*x2 + h*(-f*y1*z2 - x1*w2);
    rels[4]  := z1*y2 + h*(-x1*z2 + y1*w2);
    rels[5]  := w1*x2 + h*(-x1*z2 - f*y1*w2);
    rels[6]  := w1*y2 + h*(y1*z2 - x1*w2);

    rels[7]  := y2*x3 - x2*y3;
    rels[8]  := w2*z3 + z2*w3;
    rels[9]  := z2*x3 + h*(-f*y2*z3 - x2*w3);
    rels[10] := z2*y3 + h*(-x2*z3 + y2*w3);
    rels[11] := w2*x3 + h*(-x2*z3 - f*y2*w3);
    rels[12] := w2*y3 + h*(y2*z3 - x2*w3);

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraX" then
    # h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(Rationals, ["h"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    h := ind[1];;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 - x1*y2;
    rels[2]  := w1*z2 + z1*w2;
    rels[3]  := z1*x2 + h*(-x1*w2);
    rels[4]  := z1*y2 + h*(-x1*w2 - y1*w2);
    rels[5]  := w1*x2 + h*(-x1*z2);
    rels[6]  := w1*y2 + h*(-x1*z2 - y1*z2);

    rels[7]  := y2*x3 - x2*y3;
    rels[8]  := w2*z3 + z2*w3;
    rels[9]  := z2*x3 + h*(-x2*w3);
    rels[10] := z2*y3 + h*(-x2*w3 - y2*w3);
    rels[11] := w2*x3 + h*(-x2*z3);
    rels[12] := w2*y3 + h*(-x2*z3 - y2*z3);

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraY" then
    # f, h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(Rationals, ["f","h"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    f := ind[1];; h := ind[2];;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 - x1*y2;
    rels[2]  := w1*z2 + z1*w2;
    rels[3]  := z1*x2 + h*(-x1*z2);
    rels[4]  := z1*y2 + h*(-f*x1*z2 + y1*z2 - x1*w2);
    rels[5]  := w1*x2 + h*(-x1*w2);
    rels[6]  := w1*y2 + h*(-x1*z2 - f*x1*w2 + y1*w2);

    rels[7]  := y2*x3 - x2*y3;
    rels[8]  := w2*z3 + z2*w3;
    rels[9]  := z2*x3 + h*(-x2*z3);
    rels[10] := z2*y3 + h*(-f*x2*z3 + y2*z3 - x2*w3);
    rels[11] := w2*x3 + h*(-x2*w3);
    rels[12] := w2*y3 + h*(-x2*z3 - f*x2*w3 + y2*w3);

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

# ============================================================
if type = "AlgebraZ" then
    # f, h variable
    # x1->x, x2->y, x3->z, x4->w

    K := FunctionField(Rationals, ["f","h"]);
    kQ := PathAlgebra(K, Q);
    ind := IndeterminatesOfFunctionField(K);
    f := ind[1];; h := ind[2];;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 + x1*y2;
    rels[2]  := w1*z2 - z1*w2;
    rels[3]  := z1*x2 + h*(-x1*z2 - y1*w2);
    rels[4]  := z1*y2 + h*(-y1*z2 - x1*w2);
    rels[5]  := w1*x2 + h*(-f*y1*z2 + x1*w2);
    rels[6]  := w1*y2 + h*(-f*x1*z2 + y1*w2);

    rels[7]  := y2*x3 + x2*y3;
    rels[8]  := w2*z3 - z2*w3;
    rels[9]  := z2*x3 + h*(-x2*z3 - y2*w3);
    rels[10] := z2*y3 + h*(-y2*z3 - x2*w3);
    rels[11] := w2*x3 + h*(-f*y2*z3 + x2*w3);
    rels[12] := w2*y3 + h*(-f*x2*z3 + y2*w3);

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;
