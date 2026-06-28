type := "Sinfinity";
# type := "KirkmanS";

if type = "LBA" then # "Central Extensions of Three Dimensional Artin-Schelter Regular Algebras", by Le Bruy

    # Define coefficient field
    # K := FunctionField(Rationals, 
        #  ["a","b"]);

    K:=Rationals;
    # Path algebra
    kQ := PathAlgebra(K, Q);

    

    c := 1;

    # Extract indeterminates
    # ind := IndeterminatesOfFunctionField(K);

    a := 1/3;;
    b := 2/3;;
    alpha1:=1/3;;
    alpha2:=3/4;;
    alpha3:=1/3;;

    M11:=22/7;;
    M12:=42/13;;
    M21:=42/13;;
    M22:=7/13;;
    M31:=15/39;;
    M13:=15/39;;
    M23:=17/21;;
    M32:=17/21;;
    M33:=19/27;;


    # Arrows
    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;

    rels := [];

    rels[1] := c*x1*x2 + a*y1*z2 + b*z1*y2 
               + M11*x1*w2 + M12*y1*w2 + M13*z1*w2 + alpha1*w1*w2;

    rels[2] := c*y1*y2 + a*z1*x2 + b*x1*z2
               + M21*x1*w2 + M22*y1*w2 + M23*z1*w2 + alpha2*w1*w2;

    rels[3] := c*z1*z2 + a*x1*y2 + b*y1*x2
               + M31*x1*w2 + M32*y1*w2 + M33*z1*w2 + alpha3*w1*w2;

    rels[4] := x1*w2 - w1*x2;
    rels[5] := y1*w2 - w1*y2;
    rels[6] := z1*w2 - w1*z2;

    rels[7] := c*x2*x3 + a*y2*z3 + b*z2*y3
               + M11*x2*w3 + M12*y2*w3 + M13*z2*w3 + alpha1*w2*w3;

    rels[8] := c*y2*y3 + a*z2*x3 + b*x2*z3
               + M21*x2*w3 + M22*y2*w3 + M23*z2*w3 + alpha2*w2*w3;

    rels[9] := c*z2*z3 + a*x2*y3 + b*y2*x3
               + M31*x2*w3 + M32*y2*w3 + M33*z2*w3 + alpha3*w2*w3;

    rels[10] := x2*w3 - w2*x3;
    rels[11] := y2*w3 - w2*y3;
    rels[12] := z2*w3 - w2*z3;

    I := TwoSidedIdeal(kQ, rels);

    A := kQ/I;

fi;


if type="OreTypeA3" then

    K := FunctionField(CyclotomicField(3), ["a"]);
    kQ := PathAlgebra(K, Q);

    ind := IndeterminatesOfFunctionField(K);
    a := ind[1];; 
    b := 4/7;;

    c := 2/5;; d := 3/7;;
    p := E(3) * One(K);;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := c*x1*x2 + a*y1*z2 + b*z1*y2;
    rels[2]  := c*y1*y2 + a*z1*x2 + b*x1*z2;
    rels[3]  := c*z1*z2 + a*x1*y2 + b*y1*x2;
    rels[4]  := w1*x2 - d*z1*w2;
    rels[5]  := w1*y2 - d*p*x1*w2;
    rels[6]  := w1*z2 - d*(p^2)*y1*w2;

    rels[7]  := c*x2*x3 + a*y2*z3 + b*z2*y3;
    rels[8]  := c*y2*y3 + a*z2*x3 + b*x2*z3;
    rels[9]  := c*z2*z3 + a*x2*y3 + b*y2*x3;
    rels[10] := w2*x3 - d*z2*w3;
    rels[11] := w2*y3 - d*p*x2*w3;
    rels[12] := w2*z3 - d*(p^2)*y2*w3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

if type="KQ[x,y,z,w]" then

    K := FunctionField(Rationals, ["q14","q24"]);
    kQ := PathAlgebra(K, Q);

    ind := IndeterminatesOfFunctionField(K);
    q14 := ind[1];; q24 := ind[2];;

    q12 := 2/3;; q13 := 3/5;; q23 := 4/7;; q34 := 5/9;;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 - q12*x1*y2;
    rels[2]  := z1*x2 - q13*x1*z2;
    rels[3]  := w1*x2 - q14*x1*w2;
    rels[4]  := z1*y2 - q23*y1*z2;
    rels[5]  := w1*y2 - q24*y1*w2;
    rels[6]  := w1*z2 - q34*z1*w2;

    rels[7]  := y2*x3 - q12*x2*y3;
    rels[8]  := z2*x3 - q13*x2*z3;
    rels[9]  := w2*x3 - q14*x2*w3;
    rels[10] := z2*y3 - q23*y2*z3;
    rels[11] := w2*y3 - q24*y2*w3;
    rels[12] := w2*z3 - q34*z2*w3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

if type="LBATwist" then 

    K:=Rationals;;
    # K := FunctionField(Rationals, ["q"]);
    kQ := PathAlgebra(K, Q);

    # ind := IndeterminatesOfFunctionField(K);
    q := 1/3;; b := 1/13;;

    c0 := 1/3;; c1 := 1/2;; c2 := 2/3;;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 + q*x1*y2 + b*w1*w2 + c0*z1*z2;
    rels[2]  := w1*y2 + q*y1*w2 + b*x1*x2 + c1*z1*z2;
    rels[3]  := x1*w2 + q*w1*x2 + b*y1*y2 + c2*z1*z2;
    rels[4]  := w1*z2 + z1*w2;
    rels[5]  := x1*z2 + z1*x2;
    rels[6]  := y1*z2 + z1*y2;

    rels[7]  := y2*x3 + q*x2*y3 + b*w2*w3 + c0*z2*z3;
    rels[8]  := w2*y3 + q*y2*w3 + b*x2*x3 + c1*z2*z3;
    rels[9]  := x2*w3 + q*w2*x3 + b*y2*y3 + c2*z2*z3;
    rels[10] := w2*z3 + z2*w3;
    rels[11] := x2*z3 + z2*x3;
    rels[12] := y2*z3 + z2*y3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
    fi;



if type="L(1,1,2)" then

    K := FunctionField(Rationals, ["alpha","beta"]);
    kQ := PathAlgebra(K, Q);

    ind   := IndeterminatesOfFunctionField(K);
    alpha := ind[1];; beta := ind[2];;

    p0  := 2/3;;
    p1  := 3/5;;
    lam := alpha^2/beta;;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := x1*w2 - w1*x2;
    rels[2]  := beta*y1*w2 - (1/p0)*alpha*w1*y2;
    rels[3]  := beta*y1*x2 - p1*alpha*x1*y2;
    rels[4]  := (alpha^2/beta)*z1*y2 - (1/p0)*p1*beta*y1*z2
                - alpha*(p1-p0)*(w1*w2 + lam*w1*x2 + x1*w2)
                - (1-p0^2)*w1*w2 - (p1^2-1)*x1*x2;
    rels[5]  := (alpha^2/beta)*z1*w2 - p0*alpha*w1*z2;
    rels[6]  := (alpha^2/beta)*z1*x2 - (1/p1)*alpha*x1*z2;

    rels[7]  := x2*w3 - w2*x3;
    rels[8]  := beta*y2*w3 - (1/p0)*alpha*w2*y3;
    rels[9]  := beta*y2*x3 - p1*alpha*x2*y3;
    rels[10] := (alpha^2/beta)*z2*y3 - (1/p0)*p1*beta*y2*z3
                - alpha*(p1-p0)*(w2*w3 + lam*w2*x3 + x2*w3)
                - (1-p0^2)*w2*w3 - (p1^2-1)*x2*x3;
    rels[11] := (alpha^2/beta)*z2*w3 - p0*alpha*w2*z3;
    rels[12] := (alpha^2/beta)*z2*x3 - (1/p1)*alpha*x2*z3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

if type = "LBATwist" then # "Central Extensions TWIST of Three Dimensional Artin-Schelter Regular Algebras", by Le Bruy

    # Define coefficient field
    K := FunctionField(Rationals, 
         ["q","b"]);

    # Path algebra
    kQ := PathAlgebra(K, Q);


    # Extract indeterminates
    ind := IndeterminatesOfFunctionField(K);

    q := ind[1];;
    b := ind[2];;
    c0 := 13/75;;
    c1 := 17/19;;
    c2 := 23/27;;



    # Arrows
    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;

    rels := [];

    rels[1]:=y1*x2+q*x1*y2+b*w1*w2+c0*z1*z2;
    rels[2]:=w1*y2+q*y1*w2+b*x1*x2+c1*z1*z2;
    rels[3]:=x1*w2+q*w1*x2+b*y1*y2+c2*z1*z2;
    rels[4]:=w1*z2+z1*w2;
    rels[5]:=x1*z2+z1*x2;
    rels[6]:=y1*z2+z1*y2;
    rels[7]:=y2*x3+q*x2*y3+b*w2*w3+c0*z2*z3;
    rels[8]:=w2*y3+q*y2*w3+b*x2*x3+c1*z2*z3;
    rels[9]:=x2*w3+q*w2*x3+b*y2*y3+c2*z2*z3;
    rels[10]:=w2*z3+z2*w3;
    rels[11]:=x2*z3+z2*x3;
    rels[12]:=y2*z3+z2*y3;




    I := TwoSidedIdeal(kQ, rels);

    A := kQ/I;

fi;

if type="Clifford" then

    K := Rationals;
    kQ := PathAlgebra(K, Q);

    # a_ijl values
    a121:=1/2;; a122:=1/3;; a123:=1/4;; a124:=1/5;;
    a131:=2/3;; a132:=1/4;; a133:=3/5;; a134:=1/6;;
    a141:=1/3;; a142:=2/5;; a143:=1/7;; a144:=3/8;;
    a231:=3/4;; a232:=1/5;; a233:=2/7;; a234:=1/3;;
    a241:=1/4;; a242:=3/7;; a243:=2/9;; a244:=1/5;;
    a341:=2/5;; a342:=1/6;; a343:=3/7;; a344:=1/4;;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := x1*y2 + y1*x2 + a121*x1*x2 + a122*y1*y2 + a123*z1*z2 + a124*w1*w2;
    rels[2]  := x1*z2 + z1*x2 + a131*x1*x2 + a132*y1*y2 + a133*z1*z2 + a134*w1*w2;
    rels[3]  := x1*w2 + w1*x2 + a141*x1*x2 + a142*y1*y2 + a143*z1*z2 + a144*w1*w2;
    rels[4]  := y1*z2 + z1*y2 + a231*x1*x2 + a232*y1*y2 + a233*z1*z2 + a234*w1*w2;
    rels[5]  := y1*w2 + w1*y2 + a241*x1*x2 + a242*y1*y2 + a243*z1*z2 + a244*w1*w2;
    rels[6]  := z1*w2 + w1*z2 + a341*x1*x2 + a342*y1*y2 + a343*z1*z2 + a344*w1*w2;

    rels[7]  := x2*y3 + y2*x3 + a121*x2*x3 + a122*y2*y3 + a123*z2*z3 + a124*w2*w3;
    rels[8]  := x2*z3 + z2*x3 + a131*x2*x3 + a132*y2*y3 + a133*z2*z3 + a134*w2*w3;
    rels[9]  := x2*w3 + w2*x3 + a141*x2*x3 + a142*y2*y3 + a143*z2*z3 + a144*w2*w3;
    rels[10] := y2*z3 + z2*y3 + a231*x2*x3 + a232*y2*y3 + a233*z2*z3 + a234*w2*w3;
    rels[11] := y2*w3 + w2*y3 + a241*x2*x3 + a242*y2*y3 + a243*z2*z3 + a244*w2*w3;
    rels[12] := z2*w3 + w2*z3 + a341*x2*x3 + a342*y2*y3 + a343*z2*z3 + a344*w2*w3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;



if type="Sklyanin" then

    # K := FunctionField(Rationals, 
    #      ["beta","gamm"]);
    K:=Rationals;
    # ind := IndeterminatesOfFunctionField(K);

    # alpha := ind[1];;
    beta := 1/3;;
    gamm := 2/3;;

    
    kQ := PathAlgebra(K, Q);

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels:=[];

    rels[1]:= x1*y2 - y1*x2 +((beta+gamm)/(1+beta*gamm))*(z1*w2 + w1*z2);
    rels[2]:= x1*y2 + y1*x2 - (z1*w2 - w1*z2);
    rels[3]:= x1*z2 - z1*x2 - beta*(w1*y2 + y1*w2);
    rels[4]:= x1*z2 + z1*x2 - (w1*y2 - y1*w2);
    rels[5]:= x1*w2 - w1*x2 - gamm*(y1*z2 + z1*y2);
    rels[6]:= x1*w2 + w1*x2 - (y1*z2 - z1*y2);
    rels[7]:=  x2*y3 - y2*x3 +((beta+gamm)/(1+beta*gamm))*(z2*w3 + w2*z3);
    rels[8]:=  x2*y3 + y2*x3 - (z2*w3 - w2*z3);
    rels[9]:=  x2*z3 - z2*x3 - beta*(w2*y3 + y2*w3);
    rels[10]:= x2*z3 + z2*x3 - (w2*y3 - y2*w3);
    rels[11]:= x2*w3 - w2*x3 - gamm*(y2*z3 + z2*y3);
    rels[12]:= x2*w3 + w2*x3 - (y2*z3 - z2*y3);
    
    I := TwoSidedIdeal(kQ, rels);

    A := kQ/I;
fi;


if type="SklyaninTwist" then

    # K := FunctionField(Rationals, 
    #      ["beta","gamm"]);
    K:=Rationals;
    # ind := IndeterminatesOfFunctionField(K);

    # alpha := ind[1];;
    beta := 1/3;;
    gamm := 2/3;;

    
    kQ := PathAlgebra(K, Q);

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels:=[];

    rels[1]:= x1*y2 - y1*x2 +((beta+gamm)/(1+beta*gamm))*(z1*w2 + w1*z2);
    rels[2]:= x1*y2 + y1*x2 - (z1*w2 - w1*z2);
    rels[3]:= x1*z2 - z1*x2 - beta*(w1*y2 + y1*w2);
    rels[4]:= x1*z2 + z1*x2 - (w1*y2 - y1*w2);
    rels[5]:= x1*w2 - w1*x2 + gamm*(y1*z2 + z1*y2);
    rels[6]:= x1*w2 + w1*x2 + (y1*z2 - z1*y2);
    rels[7]:=  x2*y3 - y2*x3 +((beta+gamm)/(1+beta*gamm))*(z2*w3 + w2*z3);
    rels[8]:=  x2*y3 + y2*x3 - (z2*w3 - w2*z3);
    rels[9]:=  x2*z3 - z2*x3 - beta*(w2*y3 + y2*w3);
    rels[10]:= x2*z3 + z2*x3 - (w2*y3 - y2*w3);
    rels[11]:= x2*w3 - w2*x3 + gamm*(y2*z3 + z2*y3);
    rels[12]:= x2*w3 + w2*x3 + (y2*z3 - z2*y3);
    
    I := TwoSidedIdeal(kQ, rels);

    A := kQ/I;
fi;




if type="VancliffTw" then 

    K := FunctionField(Rationals, 
         ["alpha","beta"]);

    # Path algebra
    kQ := PathAlgebra(K, Q);


    # Extract indeterminates
    ind := IndeterminatesOfFunctionField(K);

    alpha := ind[1];;
    beta := ind[2];;
    lam := 13/27;;
    



    # Arrows
    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels:= [ ];

    rels[1]:= y1*x2 - alpha*x1*y2 ;
    rels[2]:= z1*x2 - lam*x1*z2 ;
    rels[3]:= w1*x2 - alpha*lam*x1*w2 ;
    rels[4]:= w1*z2 + alpha*z1*w2 ;
    rels[5]:= w1*y2 + lam*y1*w2 ;
    rels[6]:= z1*y2 + beta*y1*z2 - (alpha*beta - lam)*x1*w2 ;

    rels[7]:= y2*x3 - alpha*x2*y3 ;
    rels[8]:= z2*x3 - lam*x2*z3 ;
    rels[9]:= w2*x3 - alpha*lam*x2*w3 ;
    rels[10]:= w2*z3 + alpha*z2*w3 ;
    rels[11]:= w2*y3 + lam*y2*w3 ;
    rels[12]:= z2*y3 + beta*y2*z3 - (alpha*beta - lam)*x2*w3 ;

    I := TwoSidedIdeal(kQ, rels);

    A := kQ/I;

fi;

if type="KirkmanR2" then

    K:=Rationals;
    kQ:=PathAlgebra(K,Q);
    
    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;

    rels:= [ ] ;
	rels[1]:= x1*y2 + y1*x2 ;
    rels[2]:= x1*z2 + w1*y2 ;
    rels[3]:= x1*w2 - z1*y2 ;
    rels[4]:= y1*z2 - w1*x2 ;
    rels[5]:= y1*w2 + z1*x2 ;
    rels[6]:= z1*w2 + w1*z2 ;
    rels[7]:= x2*y3 + y2*x3 ;
    rels[8]:= x2*z3 + w2*y3 ;
    rels[9]:= x2*w3 - z2*y3;
    rels[10]:= y2*z3 - w2*x3 ;
    rels[11]:= y2*w3 + z2*x3 ;
    rels[12]:= z2*w3 + w2*z3 ;

    I := TwoSidedIdeal(kQ, rels);

    A := kQ/I;

fi;

if type="KirkmanS" then

    K:=Rationals;
    kQ:=PathAlgebra(K,Q);
    
    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;

    rels:= [ ] ;
	rels[1]:= x1*y2 - z1*z2 ;
    rels[2]:= x1*z2 - y1*w2 ;
    rels[3]:= x1*w2 - w1*y2 ;
    rels[4]:= y1*z2 - z1*x2 ;
    rels[5]:= z1*y2 - w1*x2 ;
    rels[6]:= y1*x2 - w1*w2 ;
    rels[7]:= x2*y3 - z2*z3 ;
    rels[8]:= x2*z3 - y2*w3 ;
    rels[9]:= x2*w3 - w2*y3 ;
    rels[10]:= y2*z3 - z2*x3 ;
    rels[11]:= z2*y3 - w2*x3 ;
    rels[12]:= y2*x3 - w2*w3 ;

    I := TwoSidedIdeal(kQ, rels);

    A := kQ/I;

fi;

if type="KirkmanT" then

    K:=Rationals;
    kQ:=PathAlgebra(K,Q);
    
    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;

    rels:= [ ] ;
	rels[1]:= x1*y2 - z1*z2 ;
    rels[2]:= x1*z2 - y1*w2 ;
    rels[3]:= x1*w2 + w1*y2 ;
    rels[4]:= y1*z2 - z1*x2 ;
    rels[5]:= z1*y2 - w1*x2 ;
    rels[6]:= y1*x2 + w1*w2 ;
    rels[7]:= x2*y3 - z2*z3 ;
    rels[8]:= x2*z3 - y2*w3 ;
    rels[9]:= x2*w3 + w2*y3 ;
    rels[10]:= y2*z3 - z2*x3 ;
    rels[11]:= z2*y3 - w2*x3 ;
    rels[12]:= y2*x3 + w2*w3 ;

    I := TwoSidedIdeal(kQ, rels);

    A := kQ/I;

fi;

if type="CentralExtS1Dash" then

    K := FunctionField(Rationals, 
         ["a","alpha"]);

    # Path algebra
    kQ := PathAlgebra(K, Q);


    # Extract indeterminates
    ind := IndeterminatesOfFunctionField(K);

    a := ind[1];;
    alpha := ind[2];;
    



    # Arrows
    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels:= [ ];

    rels[1]:=y1*z2+a*1/alpha*z1*y2;
    rels[2]:=alpha*z1*x2+a*x1*z2;
    rels[3]:=z1*z2+x1*y2+a*y1*x2;
    rels[4]:=w1*x2-x1*w2;
    rels[5]:=w1*y2-y1*w2;
    rels[6]:=w1*z2-z1*w2;
    
    rels[7]:=y2*z3+a*1/alpha*z2*y3;
    rels[8]:=alpha*z2*x3+a*x2*z3;
    rels[9]:=z2*z3+x2*y3+a*y2*x3;
    rels[10]:=w2*x3-x2*w3;
    rels[11]:=w2*y3-y2*w3;
    rels[12]:=w2*z3-z2*w3;

    I := TwoSidedIdeal(kQ, rels);

    A := kQ/I;
fi;

if type="CentralExtS1" then

    K := FunctionField(Rationals, 
         ["a","alpha"]);

    # Path algebra
    kQ := PathAlgebra(K, Q);


    # Extract indeterminates
    ind := IndeterminatesOfFunctionField(K);

    a := ind[1];;
    alpha := ind[2];;
    beta := 5/13;;
    



    # Arrows
    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels:= [ ];

    rels[1]:=y1*z2+a*beta*z1*y2;
    rels[2]:=alpha*z1*x2+a*x1*z2;
    rels[3]:=x1*y2+a*y1*x2;
    rels[4]:=w1*x2-x1*w2;
    rels[5]:=w1*y2-y1*w2;
    rels[6]:=w1*z2-z1*w2;
    
    rels[7]:=y2*z3+a*beta*z2*y3;
    rels[8]:=alpha*z2*x3+a*x2*z3;
    rels[9]:=x2*y3+a*y2*x3;
    rels[10]:=w2*x3-x2*w3;
    rels[11]:=w2*y3-y2*w3;
    rels[12]:=w2*z3-z2*w3;

    I := TwoSidedIdeal(kQ, rels);

    A := kQ/I;
fi;


if type="CentralExtS2" then

    K := FunctionField(Rationals, 
         ["alpha"]);

    # Path algebra
    kQ := PathAlgebra(K, Q);


    # Extract indeterminates
    ind := IndeterminatesOfFunctionField(K);

    alpha := ind[1];;
    



    # Arrows
    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels:= [ ];

    rels[1]:=z1*x2+1/alpha*x1*z2;
    rels[2]:=z1*y2-1/alpha*y1*z2;
    rels[3]:=x1*x2-y1*y2;
    rels[4]:=w1*x2-x1*w2;
    rels[5]:=w1*y2-y1*w2;
    rels[6]:=w1*z2-z1*w2;
    
    rels[7]:=z2*x3+1/alpha*x2*z3;
    rels[8]:=z2*y3-1/alpha*y2*z3;
    rels[9]:=x2*x3-y2*y3;
    rels[10]:=w2*x3-x2*w3;
    rels[11]:=w2*y3-y2*w3;
    rels[12]:=w2*z3-z2*w3;

    I := TwoSidedIdeal(kQ, rels);

    A := kQ/I;
fi;

if type="R3a" then 

    # K := FunctionField(Rationals, ["a"]);
    K:=Rationals;
    kQ := PathAlgebra(K, Q);

    # ind := IndeterminatesOfFunctionField(K);
    a := 0;;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 - x1*y2 + 4*a*x1*x2;
    rels[2]  := y1*z2 - z1*y2 + 4*(a+1)*y1*y2 - 8*(a+1)*(a+2)*x1*y2 - 4*(a+2)*x1*z2;
    rels[3]  := x1*z2 - z1*x2 + 4*a*x1*y2 + 8*a^2*x1*x2 - 8*a*x1*x2;

    rels[4]  := y1*w2 - w1*y2 + 4*(a+1)*y1*z2 + 8*a*(a+1)*y1*y2 - 4*(a+3)*x1*w2
                - (64/3)*a*(a+1)*(a+2)*x1*y2 - 16*(a+1)*(a+2)*x1*z2;
    rels[5]  := x1*w2 - w1*x2 + 4*a*x1*z2 - 8*(a-a^2)*x1*y2 - (64/6)*(-a^3+3*a^2-2*a)*x1*x2;
    rels[6]  := z1*w2 - w1*z2 + 4*(a+2)*z1*z2 - 8*(a+2)*(a+3)*y1*z2 - 4*(a+3)*y1*w2
                + (64/6)*(a+2)*(a+3)*(a+4)*x1*z2 + 8*(a+3)*(a+4)*x1*w2;

    rels[7]  := y2*x3 - x2*y3 + 4*a*x2*y3;
    rels[8]  := y2*z3 - z2*y3 + 4*(a+1)*y2*y3 - 8*(a+1)*(a+2)*x2*y3 - 4*(a+2)*x2*z3;
    rels[9]  := x2*z3 - z2*x3 + 4*a*x2*y3 + 8*a^2*x2*x3 - 8*a*x2*x3;

    rels[10] := y2*w3 - w2*y3 + 4*(a+1)*y2*z3 + 8*a*(a+1)*y2*y3 - 4*(a+3)*x2*w3
                - (64/3)*a*(a+1)*(a+2)*x2*y3 - 16*(a+1)*(a+2)*x2*z3;
    rels[11] := x2*w3 - w2*x3 + 4*a*x2*z3 - 8*(a-a^2)*x2*y3 - (64/6)*(-a^3+3*a^2-2*a)*x2*x3;
    rels[12] := z2*w3 - w2*z3 + 4*(a+2)*z2*z3 - 8*(a+2)*(a+3)*y2*z3 - 4*(a+3)*y2*w3
                + (64/6)*(a+2)*(a+3)*(a+4)*x2*z3 + 8*(a+3)*(a+4)*x2*w3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

if type="Sinfinity" then

    K := Rationals;
    kQ := PathAlgebra(K, Q);

    alpha := 3/2;; beta := 2;; gamm := (-alpha - beta)/(1 + alpha*beta);;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := x1*y2 - y1*x2 - alpha*(z1*w2+w1*z2);
    rels[2]  := x1*y2 + y1*x2 - z1*w2+w1*z2;
    rels[3]  := x1*z2 - z1*x2 - beta*(w1*y2+y1*w2);
    rels[4]  := x1*z2 + z1*x2 - w1*y2+y1*w2;
    rels[5]  := -x1*x2 + y1*y2 + z1*z2 + w1*w2;
    rels[6]  := y1*y2 + ((1+alpha)/(1-beta))*z1*z2 + ((1-alpha)/(1+gamm))*w1*w2;

    rels[2]  := x2*y3 - y2*x3 - alpha*(z2*w3+w2*z3);
    rels[3]  := x2*y3 + y2*x3 - z2*w3+w2*z3;
    rels[4]  := x2*z3 - z2*x3 - beta*(w2*y3+y2*w3);
    rels[5]  := x2*z3 + z2*x3 - w2*y3+y2*w3;
    rels[6]  := -x2*x3 + y2*y3 + z2*z3 + w2*w3;
    rels[7]  := y2*y3 + ((1+alpha)/(1-beta))*z2*z3 + ((1-alpha)/(1+gamm))*w2*w3;
    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;

fi;



if type="SinfinityTwist" then

    K := Rationals;
    kQ := PathAlgebra(K, Q);

    alpha := 2/3;; beta := 3/5;; gamm := (-alpha - beta)/(1 + alpha*beta);;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 - x1*y2 - alpha*z1*w2;
    rels[2]  := y1*x2 + x1*y2 - z1*w2;
    rels[3]  := z1*x2 - x1*z2 - beta*w1*y2;
    rels[4]  := z1*x2 + x1*z2 - w1*y2;
    rels[5]  := x1*x1 + y1*y1 + z1*z1 + w1*w1;
    rels[6]  := y1*y1 + ((1+alpha)/(1-beta))*z1*z1 + ((1-alpha)/(1+gamm))*w1*w1;

    rels[7]  := y2*x3 - x2*y3 - alpha*z2*w3;
    rels[8]  := y2*x3 + x2*y3 - z2*w3;
    rels[9]  := z2*x3 - x2*z3 - beta*w2*y3;
    rels[10] := z2*x3 + x2*z3 - w2*y3;
    rels[11] := x2*x2 + y2*y2 + z2*z2 + w2*w2;
    rels[12] := y2*y2 + ((1+alpha)/(1-beta))*z2*z2 + ((1-alpha)/(1+gamm))*w2*w2;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;

fi;



if type="OreExtTypeA1" then 

    K := FunctionField(CyclotomicField(3), ["a"]);
    kQ := PathAlgebra(K, Q);

    ind := IndeterminatesOfFunctionField(K);
    a := ind[1];; 
    b := 14/29*One(K); 
    c:=11/13*One(K);
    d:=17/23*One(K);
    p := E(3) * One(K);;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := c*x1*x2 + a*y1*z2 + b*z1*y2;
    rels[2]  := c*y1*y2 + a*z1*x2 + b*x1*z2;
    rels[3]  := c*z1*z2 + a*x1*y2 + b*y1*x2;
    rels[4]  := w1*x2 - d*x1*w2;
    rels[5]  := w1*y2 - d*p*y1*w2;
    rels[6]  := w1*z2 - d*p*p*z1*w2;

    rels[7]  := c*x2*x3 + a*y2*z3 + b*z2*y3;
    rels[8]  := c*y2*y3 + a*z2*x3 + b*x2*z3;
    rels[9]  := c*z2*z3 + a*x2*y3 + b*y2*x3;
    rels[10] := w2*x3 - d*x2*w3;
    rels[11] := w2*y3 - d*p*y2*w3;
    rels[12] := w2*z3 - d*p*p*z2*w3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

if type="OreExtTypeA2" then 

    K := FunctionField(CyclotomicField(3), ["a"]);
    kQ := PathAlgebra(K, Q);

    ind := IndeterminatesOfFunctionField(K);
    a := ind[1];; 
    b := 14/29*One(K);; 
    c:=11/13*One(K);
    d:=17/23*One(K);
    p := E(3) * One(K);;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := c*x1*x2 + a*y1*z2 + b*z1*y2;
    rels[2]  := c*y1*y2 + a*z1*x2 + b*x1*z2;
    rels[3]  := c*z1*z2 + a*x1*y2 + b*y1*x2;
    rels[4]  := w1*x2 - d*y1*w2;
    rels[5]  := w1*y2 - d*p*z1*w2;
    rels[6]  := w1*z2 - d*p*p*x1*w2;

    rels[7]  := c*x2*x3 + a*y2*z3 + b*z2*y3;
    rels[8]  := c*y2*y3 + a*z2*x3 + b*x2*z3;
    rels[9]  := c*z2*z3 + a*x2*y3 + b*y2*x3;
    rels[10] := w2*x3 - d*y2*w3;
    rels[11] := w2*y3 - d*p*z2*w3;
    rels[12] := w2*z3 - d*p*p*x2*w3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

if type="A5" then

    K := FunctionField(CyclotomicField(4),
         ["d"]);
    kQ := PathAlgebra(K, Q);


    ind := IndeterminatesOfFunctionField(K);

    d := ind[1];;
    a1:=23/47*One(K);
    a4:=37/59*One(K);
    a7:=11/17*One(K);
    i:=E(4)*One(K);
    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels:= [ ];
    
    rels[1]:=d*x1*w2+w1*x2;
    rels[2]:=d*y1*w2-w1*y2;
    rels[3]:=a1*x1*x2+a4*d*d*y1*y2;
    rels[4]:=d*z1*w2-i*w1*z2;
    rels[5]:=a4*x1*y2+a4*y1*x2-a7*z1*z2;    
    rels[6]:=y1*z2+z1*y2;
    rels[7]:=d*x2*w3+w2*x3;
    rels[8]:=d*y2*w3-w2*y3;
    rels[9]:=a1*x2*x3+a4*d*d*y2*y3;
    rels[10]:=d*z2*w3-i*w2*z3;
    rels[11]:=a4*x2*y3+a4*y2*x3-a7*z2*z3;    
    rels[12]:=y2*z3+z2*y3;
    
    I := TwoSidedIdeal(kQ, rels);

    A := kQ/I;
fi;


if type="CassidyVancliff1" then

    K := FunctionField(CyclotomicField(4),
         ["gamm"]);

    # Path algebra
    kQ := PathAlgebra(K, Q);


    # Extract indeterminates
    ind := IndeterminatesOfFunctionField(K);

    gamm := ind[1];;
    
    alpha := E(4) * One(K);;
    beta  := E(4) * One(K);;


    # Arrows
    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels:= [ ];

    rels[1]:=w1*x2-alpha*x1*w2;
    rels[2]:=z1*y2-beta*y1*z2;
    rels[3]:=z1*z2-x1*x2;
    rels[4]:=w1*w2-y1*y2;
    rels[5]:=z1*x2-x1*z2+y1*y2;
    rels[6]:=w1*y2-y1*w2+gamm*gamm*x1*x2;
    rels[7]:=w2*x3-alpha*x2*w3;
    rels[8]:=z2*y3-beta*y2*z3;
    rels[9]:=z2*z3-x2*x3;
    rels[10]:=w2*w3-y2*y3;
    rels[11]:=z2*x3-x2*z3+y2*y3;
    rels[12]:=w2*y3-y2*w3+gamm*gamm*x2*x3;
    
    I := TwoSidedIdeal(kQ, rels);

    A := kQ/I;



fi;

if type="CassidyVancliff2" then

    K := FunctionField(Rationals,
         ["a","c","d"]);

    # Path algebra
    kQ := PathAlgebra(K, Q);


    # Extract indeterminates
    ind := IndeterminatesOfFunctionField(K);

    a := ind[1];;
    c:=ind[2];
    d:=ind[3];
    b := 0 * One(K);;


    # Arrows
    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels:= [ ];

    rels[1]:=z1*x2+x1*z2-d*y1*y2;
    rels[2]:=w1*x2+x1*w2-b*z1*z2;
    rels[3]:=y1*z2-z1*y2;
    rels[4]:=w1*w2-y1*y2;
    rels[5]:=w1*y2+y1*w2-z1*z2;
    rels[6]:=a*z1*z2+c*y1*y2-x1*x2;
    
    rels[7]:=z2*x3+x2*z3-d*y2*y3;
    rels[8]:=w2*x3+x2*w3-b*z2*z3;
    rels[9]:=y2*z3-z2*y3;
    rels[10]:=w2*w3-y2*y3;
    rels[11]:=w2*y3+y2*w3-z2*z3;
    rels[12]:=a*z2*z3+c*y2*y3-x2*x3;
    
    I := TwoSidedIdeal(kQ, rels);

    A := kQ/I;
fi;

if type="CassidyVancliff3" then

    K := CyclotomicField(4);

    # Path algebra
    kQ := PathAlgebra(K, Q);


    # Extract indeterminates
    a:=1;
    b:=1;
    c:=1;
    d:=1;


    # Arrows
    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels:= [ ];

    rels[1]:=x1*z2+a*z1*x2;
    rels[2]:=x1*w2+b*w1*x2;
    rels[3]:=z1*w2+d*w1*z2;
    rels[4]:=w1*w2-y1*y2;
    rels[5]:=y1*z2+z1*y2+w1*w2;
    rels[6]:=y1*w2+c*w1*y2+x1*x2;
    
    rels[7]:=x2*z3+a*z2*x3;
    rels[8]:=x2*w3+b*w2*x3;
    rels[9]:=z2*w3+d*w2*z3;
    rels[10]:=w2*w3-y2*y3;
    rels[11]:=y2*z3+z2*y3+w2*w3;
    rels[12]:=y2*w3+c*w2*y3+x2*x3;
    I := TwoSidedIdeal(kQ, rels);

    A := kQ/I;



fi;

if type="F1" then #F(0,-1,-1,2)(-1,0,-1,2)(0,0,0,0)


    K := FunctionField(Rationals,
         ["a","b"]);

    # Path algebra
    kQ := PathAlgebra(K, Q);


    # Extract indeterminates
    ind := IndeterminatesOfFunctionField(K);

    a := ind[1];;
    b:=ind[2];


    # Arrows
    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels:= [ ];

    rels[1]:=y1*x2-a*x1*y2;
    rels[2]:=z1*x2-(b*b/a)*x1*z2-w1*w2;
    rels[3]:=w1*x2-b*x1*w2;
    rels[4]:=z1*y2-b*b*a*y1*z2-w1*w2;
    rels[5]:=w1*y2-b*y1*w2;
    rels[6]:=w1*z2-(1/b)*z1*w2;
    rels[7]:=y2*x3-a*x2*y3;
    rels[8]:=z2*x3-(b*b/a)*x2*z3-w2*w3;
    rels[9]:=w2*x3-b*x2*w3;
    rels[10]:=z2*y3-b*b*a*y2*z3-w2*w3;
    rels[11]:=w2*y3-b*y2*w3;
    rels[12]:=w2*z3-(1/b)*z2*w3;


    I := TwoSidedIdeal(kQ, rels);

    A := kQ/I;

fi;

if type="F2" then 

    K := FunctionField(Rationals, ["q12","q24"]);
    kQ := PathAlgebra(K, Q);

    ind := IndeterminatesOfFunctionField(K);
    q12 := ind[1];; q24 := ind[2];; s := 1;;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 - q12*x1*y2 - z1*z2;
    rels[2]  := z1*x2 - q24^(-6)*q12*x1*z2;
    rels[3]  := w1*x2 - q24^(-3)*x1*w2;
    rels[4]  := z1*y2 - q12*q24^(-6)*y1*z2 - s*w1*w2;
    rels[5]  := w1*y2 - q24*y1*w2;
    rels[6]  := w1*z2 - q24^(-1)*z1*w2;

    rels[7]  := y2*x3 - q12*x2*y3 - z2*z3;
    rels[8]  := z2*x3 - q24^(-6)*q12*x2*z3;
    rels[9]  := w2*x3 - q24^(-3)*x2*w3;
    rels[10] := z2*y3 - q12*q24^(-6)*y2*z3 - s*w2*w3;
    rels[11] := w2*y3 - q24*y2*w3;
    rels[12] := w2*z3 - q24^(-1)*z2*w3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;

fi;

if type="Casines1" then

    K := FunctionField(Rationals, ["a"]);

    kQ := PathAlgebra(K, Q);

    ind := IndeterminatesOfFunctionField(K);
    a := ind[1];; b := 1;; 
    c := 11/17;; d := 13/19;;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := w1*z2 - z1*w2 - a*x1*y2;
    rels[2]  := w1*y2 - b*z1*z2 + y1*w2;
    rels[3]  := w1*x2 - c*z1*z2 + x1*w2;
    rels[4]  := z1*y2 - y1*z2 + (b*d/c)*y1*w2 - (b*b*d/(c*c))*x1*w2;
    rels[5]  := z1*x2 - x1*z2 + (b*d/c)*x1*w2 - d*y1*w2;
    rels[6]  := y1*x2 + x1*y2;

    rels[7]  := w2*z3 - z2*w3 - a*x2*y3;
    rels[8]  := w2*y3 - b*z2*z3 + y2*w3;
    rels[9]  := w2*x3 - c*z2*z3 + x2*w3;
    rels[10] := z2*y3 - y2*z3 + (b*d/c)*y2*w3 - (b*b*d/(c*c))*x2*w3;
    rels[11] := z2*x3 - x2*z3 + (b*d/c)*x2*w3 - d*y2*w3;
    rels[12] := y2*x3 + x2*y3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;


################
if type="Lie11" then

    K := Rationals;
    kQ := PathAlgebra(K, Q);

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*z2 - z1*y2 - x1*w2;
    rels[2]  := x1*y2 - y1*x2;
    rels[3]  := x1*z2 - z1*x2;
    rels[4]  := x1*w2 - w1*x2;
    rels[5]  := y1*w2 - w1*y2;
    rels[6]  := z1*w2 - w1*z2;

    rels[7]  := y2*z3 - z2*y3 - x2*w3;
    rels[8]  := x2*y3 - y2*x3;
    rels[9]  := x2*z3 - z2*x3;
    rels[10] := x2*w3 - w2*x3;
    rels[11] := y2*w3 - w2*y3;
    rels[12] := z2*w3 - w2*z3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

if type="Lie2" then 

    K := Rationals;
    kQ := PathAlgebra(K, Q);

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := x1*y2 - y1*x2 - x1*w2;
    rels[2]  := x1*z2 - z1*x2;
    rels[3]  := y1*z2 - z1*y2;
    rels[4]  := x1*w2 - w1*x2;
    rels[5]  := y1*w2 - w1*y2;
    rels[6]  := z1*w2 - w1*z2;

    rels[7]  := x2*y3 - y2*x3 - x2*w3;
    rels[8]  := x2*z3 - z2*x3;
    rels[9]  := y2*z3 - z2*y3;
    rels[10] := x2*w3 - w2*x3;
    rels[11] := y2*w3 - w2*y3;
    rels[12] := z2*w3 - w2*z3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;

fi;

if type="SheltonTingy" then

    K := CyclotomicField(4);
    kQ := PathAlgebra(K, Q);

    p := E(4);;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := z1*x2 - x1*z2 + y1*y2;
    rels[2]  := p*w1*x2 + x1*w2;
    rels[3]  := w1*y2 - y1*w2 + z1*z2;
    rels[4]  := p*z1*y2 + y1*z2;
    rels[5]  := x1*x2 - z1*z2;
    rels[6]  := y1*y2 - w1*w2;

    rels[7]  := z2*x3 - x2*z3 + y2*y3;
    rels[8]  := p*w2*x3 + x2*w3;
    rels[9]  := w2*y3 - y2*w3 + z2*z3;
    rels[10] := p*z2*y3 + y2*z3;
    rels[11] := x2*x3 - z2*z3;
    rels[12] := y2*y3 - w2*w3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

if type="Lie3" then

    K := FunctionField(Rationals, ["alpha"]);
    kQ := PathAlgebra(K, Q);

    ind   := IndeterminatesOfFunctionField(K);
    alpha := ind[1];;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := x1*y2 - y1*x2;
    rels[2]  := x1*z2 - z1*x2 - x1*w2;
    rels[3]  := y1*z2 - z1*y2 - alpha*y1*w2;
    rels[4]  := x1*w2 - w1*x2;
    rels[5]  := y1*w2 - w1*y2;
    rels[6]  := z1*w2 - w1*z2;

    rels[7]  := x2*y3 - y2*x3;
    rels[8]  := x2*z3 - z2*x3 - x2*w3;
    rels[9]  := y2*z3 - z2*y3 - alpha*y2*w3;
    rels[10] := x2*w3 - w2*x3;
    rels[11] := y2*w3 - w2*y3;
    rels[12] := z2*w3 - w2*z3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

if type="S(2,3)" then

    K := Rationals;
    kQ := PathAlgebra(K, Q);

    alpha := 1;;

    c1:=1/2;; c2:=2/3;; c3:=1/3;;
    d1:=3/4;; d2:=1/4;; d3:=2/5;;
    beta1:=1/3;; beta2:=2/5;; beta3:=3/7;;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := alpha*(w1*x2 - x1*w2) + (beta1-alpha)*x1*x2 + beta2*y1*x2 + beta3*z1*x2
                + alpha*(c3+2)*x1*y2 - alpha*c1*x1*z2 - d1*alpha*y1*z2;
    rels[2]  := alpha*(w1*y2 - y1*w2) + beta1*x1*y2 + (beta2-alpha)*y1*y2 + beta3*z1*y2
                + alpha*(c1+2)*y1*z2 - alpha*c2*y1*x2 - d2*alpha*z1*x2;
    rels[3]  := alpha*(w1*z2 - z1*w2) + beta1*x1*z2 + beta2*y1*z2 + (beta3-alpha)*z1*z2
                + alpha*(c2+2)*z1*x2 - alpha*c3*z1*y2 - d3*alpha*x1*y2;
    rels[4]  := y1*z2 - z1*y2;
    rels[5]  := z1*x2 - x1*z2;
    rels[6]  := x1*y2 - y1*x2;

    rels[7]  := alpha*(w2*x3 - x2*w3) + (beta1-alpha)*x2*x3 + beta2*y2*x3 + beta3*z2*x3
                + alpha*(c3+2)*x2*y3 - alpha*c1*x2*z3 - d1*alpha*y2*z3;
    rels[8]  := alpha*(w2*y3 - y2*w3) + beta1*x2*y3 + (beta2-alpha)*y2*y3 + beta3*z2*y3
                + alpha*(c1+2)*y2*z3 - alpha*c2*y2*x3 - d2*alpha*z2*x3;
    rels[9]  := alpha*(w2*z3 - z2*w3) + beta1*x2*z3 + beta2*y2*z3 + (beta3-alpha)*z2*z3
                + alpha*(c2+2)*z2*x3 - alpha*c3*z2*y3 - d3*alpha*x2*y3;
    rels[10] := y2*z3 - z2*y3;
    rels[11] := z2*x3 - x2*z3;
    rels[12] := x2*y3 - y2*x3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

if type="Lie4" then

    K := FunctionField(Rationals, ["beta"]);
    kQ := PathAlgebra(K, Q);

    ind  := IndeterminatesOfFunctionField(K);
    beta := ind[1];;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := x1*y2 - y1*x2;
    rels[2]  := x1*z2 - z1*x2 - x1*w2 - beta*y1*w2;
    rels[3]  := y1*z2 - z1*y2 - y1*w2;
    rels[4]  := x1*w2 - w1*x2;
    rels[5]  := y1*w2 - w1*y2;
    rels[6]  := z1*w2 - w1*z2;

    rels[7]  := x2*y3 - y2*x3;
    rels[8]  := x2*z3 - z2*x3 - x2*w3 - beta*y2*w3;
    rels[9]  := y2*z3 - z2*y3 - y2*w3;
    rels[10] := x2*w3 - w2*x3;
    rels[11] := y2*w3 - w2*y3;
    rels[12] := z2*w3 - w2*z3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;

fi;

if type="sl2" then

    K := Rationals;
    kQ := PathAlgebra(K, Q);

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := x1*z2 - z1*x2 - 2*x1*w2;
    rels[2]  := y1*z2 - z1*y2 + 2*y1*w2;
    rels[3]  := x1*y2 - y1*x2 - z1*w2;
    rels[4]  := x1*w2 - w1*x2;
    rels[5]  := y1*w2 - w1*y2;
    rels[6]  := z1*w2 - w1*z2;

    rels[7]  := x2*z3 - z2*x3 - 2*x2*w3;
    rels[8]  := y2*z3 - z2*y3 + 2*y2*w3;
    rels[9]  := x2*y3 - y2*x3 - z2*w3;
    rels[10] := x2*w3 - w2*x3;
    rels[11] := y2*w3 - w2*y3;
    rels[12] := z2*w3 - w2*z3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

if type="CentralExtS2" then

    K := FunctionField(Rationals, ["alpha"]);
    kQ := PathAlgebra(K, Q);

    ind   := IndeterminatesOfFunctionField(K);
    alpha := ind[1];;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := z1*x2 + (1/alpha)*x1*z2;
    rels[2]  := z1*y2 - (1/alpha)*y1*z2;
    rels[3]  := x1*x2 - y1*y2;
    rels[4]  := w1*x2 - x1*w2;
    rels[5]  := w1*y2 - y1*w2;
    rels[6]  := w1*z2 - z1*w2;

    rels[7]  := z2*x3 + (1/alpha)*x2*z3;
    rels[8]  := z2*y3 - (1/alpha)*y2*z3;
    rels[9]  := x2*x3 - y2*y3;
    rels[10] := w2*x3 - x2*w3;
    rels[11] := w2*y3 - y2*w3;
    rels[12] := w2*z3 - z2*w3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

if type="OreTypeB1" then

    K := FunctionField(Rationals, ["a","d"]);
    kQ := PathAlgebra(K, Q);

    ind := IndeterminatesOfFunctionField(K);
    a := ind[1];; d := ind[2];;
    p := 1;;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := x1*y2 + y1*x2 + y1*y2 - z1*z2;
    rels[2]  := x1*x2 + y1*x2 + x1*y2 - a*z1*z2;
    rels[3]  := z1*x2 - x1*z2 + a*z1*y2 - a*y1*z2;
    rels[4]  := w1*x2 - d*x1*w2;
    rels[5]  := w1*y2 - d*y1*w2;
    rels[6]  := w1*z2 - p*d*z1*w2;

    rels[7]  := x2*y3 + y2*x3 + y2*y3 - z2*z3;
    rels[8]  := x2*x3 + y2*x3 + x2*y3 - a*z2*z3;
    rels[9]  := z2*x3 - x2*z3 + a*z2*y3 - a*y2*z3;
    rels[10] := w2*x3 - d*x2*w3;
    rels[11] := w2*y3 - d*y2*w3;
    rels[12] := w2*z3 - p*d*z2*w3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;

fi;

if type="OreTypeE1" then

    K := FunctionField(CyclotomicField(9), ["d"]);
    kQ := PathAlgebra(K, Q);

    ind := IndeterminatesOfFunctionField(K);
    d   := ind[1];;
    zeta := E(9) * One(K);;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := z1*x2 + zeta^8*x1*z2 + zeta^4*y1*y2;
    rels[2]  := x1*y2 + zeta^5*y1*x2 + zeta^7*z1*z2;
    rels[3]  := zeta*x1*x2 + y1*z2 + zeta^2*z1*y2;
    rels[4]  := w1*x2 - d*x1*w2;
    rels[5]  := w1*y2 - d*y1*w2;
    rels[6]  := w1*z2 - d*z1*w2;

    rels[7]  := z2*x3 + zeta^8*x2*z3 + zeta^4*y2*y3;
    rels[8]  := x2*y3 + zeta^5*y2*x3 + zeta^7*z2*z3;
    rels[9]  := zeta*x2*x3 + y2*z3 + zeta^2*z2*y3;
    rels[10] := w2*x3 - d*x2*w3;
    rels[11] := w2*y3 - d*y2*w3;
    rels[12] := w2*z3 - d*z2*w3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;


fi;

if type="OreTypeE2" then

    K := FunctionField(CyclotomicField(9), ["d"]);
    kQ := PathAlgebra(K, Q);

    ind := IndeterminatesOfFunctionField(K);
    d   := ind[1];;
    zeta := E(9) * One(K);;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := z1*x2 + zeta^8*x1*z2 + zeta^4*y1*y2;
    rels[2]  := x1*y2 + zeta^5*y1*x2 + zeta^7*z1*z2;
    rels[3]  := zeta*x1*x2 + y1*z2 + zeta^2*z1*y2;
    rels[4]  := w1*x2 - d*x1*w2;
    rels[5]  := w1*y2 - d*zeta^3*y1*w2;
    rels[6]  := w1*z2 - zeta^6*d*z1*w2;

    rels[7]  := z2*x3 + zeta^8*x2*z3 + zeta^4*y2*y3;
    rels[8]  := x2*y3 + zeta^5*y2*x3 + zeta^7*z2*z3;
    rels[9]  := zeta*x2*x3 + y2*z3 + zeta^2*z2*y3;
    rels[10] := w2*x3 - d*x2*w3;
    rels[11] := w2*y3 - d*zeta^3*y2*w3;
    rels[12] := w2*z3 - d*zeta^6*z2*w3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;

fi;

if type="OreTypeH1" then 

    K := FunctionField(CyclotomicField(4), ["d"]);
    kQ := PathAlgebra(K, Q);

    ind := IndeterminatesOfFunctionField(K);
    d    := ind[1];;
    zeta := E(4) * One(K);;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := x1*x2 - y1*y2;
    rels[2]  := x1*y2 - y1*x2 + zeta*z1*z2;
    rels[3]  := y1*z2 - zeta*z1*y2;
    rels[4]  := w1*x2 - d*x1*w2;
    rels[5]  := w1*y2 + d*y1*w2;
    rels[6]  := w1*z2 - zeta*d*z1*w2;

    rels[7]  := x2*x3 - y2*y3;
    rels[8]  := x2*y3 - y2*x3 + zeta*z2*z3;
    rels[9]  := y2*z3 - zeta*z2*y3;
    rels[10] := w2*x3 - d*x2*w3;
    rels[11] := w2*y3 + d*y2*w3;
    rels[12] := w2*z3 - zeta*d*z2*w3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;

fi;

if type="OreTypeH2" then 

    K := FunctionField(CyclotomicField(4), ["d"]);
    kQ := PathAlgebra(K, Q);

    ind := IndeterminatesOfFunctionField(K);
    d    := ind[1];;
    zeta := E(4) * One(K);;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := []; 

    rels[1]  := x1*x2 - y1*y2;
    rels[2]  := x1*y2 - y1*x2 + zeta*z1*z2;
    rels[3]  := y1*z2 - zeta*z1*y2;
    rels[4]  := w1*x2 - d*x1*w2;
    rels[5]  := w1*y2 - d*y1*w2;
    rels[6]  := w1*z2 + d*z1*w2;

    rels[7]  := x2*x3 - y2*y3;
    rels[8]  := x2*y3 - y2*x3 + zeta*z2*z3;
    rels[9]  := y2*z3 - zeta*z2*y3;
    rels[10] := w2*x3 - d*x2*w3;
    rels[11] := w2*y3 - d*y2*w3;
    rels[12] := w2*z3 + d*z2*w3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;

fi;

if type="CentralExtH" then

    K := CyclotomicField(4);
    kQ := PathAlgebra(K, Q);

    zeta := E(4);;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := x1*x2 - y1*y2 + x1*w2;
    rels[2]  := x1*y2 - y1*x2 + zeta*z1*z2;
    rels[3]  := y1*z2 - zeta*z1*y2;
    rels[4]  := w1*x2 - x1*w2;
    rels[5]  := w1*y2 - y1*w2;
    rels[6]  := w1*z2 - z1*w2;

    rels[7]  := x2*x3 - y2*y3 + x2*w3;
    rels[8]  := x2*y3 - y2*x3 + zeta*z2*z3;
    rels[9]  := y2*z3 - zeta*z2*y3;
    rels[10] := w2*x3 - x2*w3;
    rels[11] := w2*y3 - y2*w3;
    rels[12] := w2*z3 - z2*w3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

if type="OreTypeS2" then 

    K := FunctionField(Rationals, ["a","alpha"]);
    kQ := PathAlgebra(K, Q);

    ind   := IndeterminatesOfFunctionField(K);
    a     := ind[1];; alpha := ind[2];;

    p := 1;;
    d := 17/37;;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := z1*x2 + (1/alpha)*x1*z2;
    rels[2]  := z1*y2 - (1/alpha)*y1*z2;
    rels[3]  := x1*x2 - y1*y2;
    rels[4]  := w1*x2 - a*d*x1*w2;
    rels[5]  := w1*y2 - p*a*d*y1*w2;
    rels[6]  := w1*z2 - d*z1*w2;

    rels[7]  := z2*x3 + (1/alpha)*x2*z3;
    rels[8]  := z2*y3 - (1/alpha)*y2*z3;
    rels[9]  := x2*x3 - y2*y3;
    rels[10] := w2*x3 - a*d*x2*w3;
    rels[11] := w2*y3 - p*a*d*y2*w3;
    rels[12] := w2*z3 - d*z2*w3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

if type="R(3,a)" then

    K := Rationals;
    kQ := PathAlgebra(K, Q);

    # ind := IndeterminatesOfFunctionField(K);
    a   := 1/3;;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := x1*y2 - y1*x2 + 4*a*x1*x2;
    rels[2]  := y1*z2 - z1*y2 + 4*(a+1)*y1*y2 - 8*(a+1)*(a+2)*x1*y2 - 4*(a+2)*x1*z2;
    rels[3]  := x1*z2 - z1*x2 + 4*a*x1*y2 + 8*a^2*x1*x2 - 8*a*x1*x2;

    rels[4]  := y1*w2 - w1*y2 + 4*(a+1)*y1*z2 + 8*a*(a+1)*y1*y2 - 4*(a+3)*x1*w2
                - (64/3)*a*(a+1)*(a+2)*x1*y2 - 16*(a+1)*(a+2)*x1*z2;
    rels[5]  := x1*w2 - w1*x2 + 4*a*x1*z2 - 8*(a-a^2)*x1*y2 - (64/6)*(-a^3+3*a^2-2*a)*x1*x2;
    rels[6]  := z1*w2 - w1*z2 + 4*(a+2)*z1*z2 - 8*(a+2)*(a+3)*y1*z2 - 4*(a+3)*y1*w2
                + (64/6)*(a+2)*(a+3)*(a+4)*x1*z2 + 8*(a+3)*(a+4)*x1*w2;

    rels[7]  := x2*y3 - y2*x3 + 4*a*x2*x3;
    rels[8]  := y2*z3 - z2*y3 + 4*(a+1)*y2*y3 - 8*(a+1)*(a+2)*x2*y3 - 4*(a+2)*x2*z3;
    rels[9]  := x2*z3 - z2*x3 + 4*a*x2*y3 + 8*a^2*x2*x3 - 8*a*x2*x3;

    rels[10] := y2*w3 - w2*y3 + 4*(a+1)*y2*z3 + 8*a*(a+1)*y2*y3 - 4*(a+3)*x2*w3
                - (64/3)*a*(a+1)*(a+2)*x2*y3 - 16*(a+1)*(a+2)*x2*z3;
    rels[11] := x2*w3 - w2*x3 + 4*a*x2*z3 - 8*(a-a^2)*x2*y3 - (64/6)*(-a^3+3*a^2-2*a)*x2*x3;
    rels[12] := z2*w3 - w2*z3 + 4*(a+2)*z2*z3 - 8*(a+2)*(a+3)*y2*z3 - 4*(a+3)*y2*w3
                + (64/6)*(a+2)*(a+3)*(a+4)*x2*z3 + 8*(a+3)*(a+4)*x2*w3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;
if type="OreTypeS1dash" then

    K := FunctionField(Rationals, ["a","alpha"]);
    kQ := PathAlgebra(K, Q);

    ind   := IndeterminatesOfFunctionField(K);
    a     := ind[1];; alpha := ind[2];;

    b := 1/3;; c := 11/17;; d := 13/23;;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*z2 + a*(1/alpha)*z1*y2;
    rels[2]  := alpha*z1*x2 + a*x1*z2;
    rels[3]  := z1*z2 + x1*y2 + a*y1*x2;
    rels[4]  := w1*x2 - b*b*d*x1*w2;
    rels[5]  := w1*y2 - c*c*d*y1*w2;
    rels[6]  := w1*z2 - b*c*d*z1*w2;

    rels[7]  := y2*z3 + a*(1/alpha)*z2*y3;
    rels[8]  := alpha*z2*x3 + a*x2*z3;
    rels[9]  := z2*z3 + x2*y3 + a*y2*x3;
    rels[10] := w2*x3 - b*b*d*x2*w3;
    rels[11] := w2*y3 - c*c*d*y2*w3;
    rels[12] := w2*z3 - b*c*d*z2*w3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

if type="CentralExtB" then

    K := Rationals;
    kQ := PathAlgebra(K, Q);

    # ind := IndeterminatesOfFunctionField(K);
    a   := 1;;

    l11 := 2/3;; l12 := 1/2;; l22 := 3/4;;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := x1*y2 + y1*x2 + y1*y2 - z1*z2 + l11*x1*w2 + l12*y1*w2;
    rels[2]  := x1*x2 + y1*x2 + x1*y2 - a*z1*z2 + l12*x1*w2 + l22*y1*w2;
    rels[3]  := z1*x2 - x1*z2 + a*z1*y2 - a*y1*z2;
    rels[4]  := w1*x2 - x1*w2;
    rels[5]  := w1*y2 - y1*w2;
    rels[6]  := w1*z2 - z1*w2;

    rels[7]  := x2*y3 + y2*x3 + y2*y3 - z2*z3 + l11*x2*w3 + l12*y2*w3;
    rels[8]  := x2*x3 + y2*x3 + x2*y3 - a*z2*z3 + l12*x2*w3 + l22*y2*w3;
    rels[9]  := z2*x3 - x2*z3 + a*z2*y3 - a*y2*z3;
    rels[10] := w2*x3 - x2*w3;
    rels[11] := w2*y3 - y2*w3;
    rels[12] := w2*z3 - z2*w3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

if type="F(-1,-1,-1,1,1)(0,0,0,0)" then 

        K := FunctionField(Rationals, ["q12","q13"]);
        kQ := PathAlgebra(K, Q);

        ind := IndeterminatesOfFunctionField(K);
        q12 := ind[1];; q13 := ind[2];; q14 := 1/3;; q23 := 2/3;;

        w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
        w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
        w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
        rels := [];

        rels[1]  := y1*x2 - q12*x1*y2 - z1*w2;
        rels[2]  := z1*x2 - q13*x1*z2;
        rels[3]  := w1*x2 - q14*x1*w2;
        rels[4]  := z1*y2 - q23*y1*z2;
        rels[5]  := w1*y2 - (1/q13)*(1/q14)*(1/q23)*y1*w2;
        rels[6]  := w1*z2 - (1/q13)*(1/q23)*z1*w2;

        rels[7]  := y2*x3 - q12*x2*y3 - z2*w3;
        rels[8]  := z2*x3 - q13*x2*z3;
        rels[9]  := w2*x3 - q14*x2*w3;
        rels[10] := z2*y3 - q23*y2*z3;
        rels[11] := w2*y3 - (1/q13)*(1/q14)*(1/q23)*y2*w3;
        rels[12] := w2*z3 - (1/q13)*(1/q23)*z2*w3;

        I := TwoSidedIdeal(kQ, rels);
        A := kQ/I;
fi;


if type="F(0,-1,-1,2)(0,0,0,0)(2,-1,-1,0)" then 


        K := FunctionField(Rationals, ["q13","q23","q34"]);
        kQ := PathAlgebra(K, Q);

        ind := IndeterminatesOfFunctionField(K);
        q13 := ind[1];; q23 := ind[2];; q34 := ind[3];;

        w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
        w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
        w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
        rels := [];

        rels[1]  := y1*x2 - (1/q13)*x1*y2;
        rels[2]  := z1*x2 - q13*x1*z2;
        rels[3]  := w1*x2 + x1*w2;
        rels[4]  := z1*y2 - q23*y1*z2 - w1*w2 - x1*x2;
        rels[5]  := w1*y2 - (1/q34)*y1*w2;
        rels[6]  := w1*z2 - q34*z1*w2;

        rels[7]  := y2*x3 - (1/q13)*x2*y3;
        rels[8]  := z2*x3 - q13*x2*z3;
        rels[9]  := w2*x3 + x2*w3;
        rels[10] := z2*y3 - q23*y2*z3 - w2*w3 - x2*x3;
        rels[11] := w2*y3 - (1/q34)*y2*w3;
        rels[12] := w2*z3 - q34*z2*w3;

        I := TwoSidedIdeal(kQ, rels);
        A := kQ/I;
fi;

if type="F(0,-1,-1,2)(-1,0,-1,2)(0,0,0,0)" then 

    K := FunctionField(Rationals, ["q13"]);
    kQ := PathAlgebra(K, Q);

    ind := IndeterminatesOfFunctionField(K);
    q13 := ind[1];;

    q23 := 2/3;; q34 := 3/5;;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 - (1/q13)*x1*y2;
    rels[2]  := z1*x2 - q13*x1*z2;
    rels[3]  := w1*x2 + x1*w2;             # using + for ±
    rels[4]  := z1*y2 - q23*y1*z2 - w1*w2 - x1*x2;
    rels[5]  := w1*y2 - (1/q34)*y1*w2;
    rels[6]  := w1*z2 - q34*z1*w2;

    rels[7]  := y2*x3 - (1/q13)*x2*y3;
    rels[8]  := z2*x3 - q13*x2*z3;
    rels[9]  := w2*x3 + x2*w3;             # using + for ±
    rels[10] := z2*y3 - q23*y2*z3 - w2*w3 - x2*x3;
    rels[11] := w2*y3 - (1/q34)*y2*w3;
    rels[12] := w2*z3 - q34*z2*w3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

if type="F(0,-1,-1,2)(0,0,0,0)" then 

    K := FunctionField(Rationals, ["q13","q14"]);
    kQ := PathAlgebra(K, Q);

    ind := IndeterminatesOfFunctionField(K);
    q13 := ind[1];; q14 := ind[2];;

    q23 := 2/3;; q24 := 3/5;;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 - q14^2*(1/q13)*x1*y2;
    rels[2]  := z1*x2 - q13*x1*z2;
    rels[3]  := w1*x2 - q14*x1*w2;
    rels[4]  := z1*y2 - q23*y1*z2 - w1*w2;
    rels[5]  := w1*y2 - q24*y1*w2;
    rels[6]  := w1*z2 - (1/q24)*z1*w2;

    rels[7]  := y2*x3 - q14^2*(1/q13)*x2*y3;
    rels[8]  := z2*x3 - q13*x2*z3;
    rels[9]  := w2*x3 - q14*x2*w3;
    rels[10] := z2*y3 - q23*y2*z3 - w2*w3;
    rels[11] := w2*y3 - q24*y2*w3;
    rels[12] := w2*z3 - (1/q24)*z2*w3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;

if type="CentralExtTypeE" then 

        K := CyclotomicField(9);
        kQ := PathAlgebra(K, Q);

        zeta := E(9);

        w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
        w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
        w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
        rels := [];

        rels[1]  := z1*x2 + (zeta^8)*x1*z2 + (zeta^4)*y1*y2;
        rels[2]  := x1*y2 + (zeta^5)*y1*x2 + (zeta^7)*z1*z2;
        rels[3]  := zeta*x1*x2 + y1*z2 + (zeta^2)*z1*y2;
        rels[4]  := w1*x2 - x1*w2;
        rels[5]  := w1*y2 - y1*w2;
        rels[6]  := w1*z2 - z1*w2;

        rels[7]  := z2*x3 + (zeta^8)*x2*z3 + (zeta^4)*y2*y3;
        rels[8]  := x2*y3 + (zeta^5)*y2*x3 + (zeta^7)*z2*z3;
        rels[9]  := zeta*x2*x3 + y2*z3 + (zeta^2)*z2*y3;
        rels[10] := w2*x3 - x2*w3;
        rels[11] := w2*y3 - y2*w3;
        rels[12] := w2*z3 - z2*w3;

        I := TwoSidedIdeal(kQ, rels);
        A := kQ/I;
fi;


if type="VancliffTwist" then 

    K := FunctionField(Rationals, ["alpha","beta","lambda"]);
    # K:=Rationals;
    kQ := PathAlgebra(K, Q);

    ind := IndeterminatesOfFunctionField(K);
    alpha := ind[1];; beta := ind[2];; lambda := ind[3]    ;;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];

    rels[1]  := y1*x2 - alpha*x1*y2;
    rels[2]  := z1*x2 - lambda*x1*z2;
    rels[3]  := w1*x2 - alpha*lambda*x1*w2;
    rels[4]  := w1*z2 + alpha*z1*w2;
    rels[5]  := w1*y2 + lambda*y1*w2;
    rels[6]  := z1*y2 + beta*y1*z2 - (alpha*beta - lambda)*x1*w2;

    rels[7]  := y2*x3 - alpha*x2*y3;
    rels[8]  := z2*x3 - lambda*x2*z3;
    rels[9]  := w2*x3 - alpha*lambda*x2*w3;
    rels[10] := w2*z3 + alpha*z2*w3;
    rels[11] := w2*y3 + lambda*y2*w3;
    rels[12] := z2*y3 + beta*y2*z3 - (alpha*beta - lambda)*x2*w3;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;


if type="SdiTwist" then 

    # K := FunctionField(Rationals, ["alpha","beta"]);
    K:=Rationals;
    alpha:=11/23;;
    beta:=13/27;;
    kQ := PathAlgebra(K, Q);

    # ind := IndeterminatesOfFunctionField(K);
    # alpha := ind[1];; beta := ind[2];;

    gamm := (-alpha - beta)/(1 + alpha*beta);;

    w1 := kQ.w1;; x1 := kQ.x1;; y1 := kQ.y1;; z1 := kQ.z1;;
    w2 := kQ.w2;; x2 := kQ.x2;; y2 := kQ.y2;; z2 := kQ.z2;;
    w3 := kQ.w3;; x3 := kQ.x3;; y3 := kQ.y3;; z3 := kQ.z3;;
    rels := [];


    rels[1]  := x1*y2 - y1*x2 - alpha*z1*w2 + alpha*w1*z2;

    rels[2]  := x1*z2 - z1*x2 - beta*w1*y2 + beta*y1*w2;

    rels[3]  := x1*z2 - z1*x2 + w1*y2 - y1*w2;

    rels[4]  := x1*w2 - w1*x2 + gamm*y1*z2 - gamm*z1*y2;

    rels[5]  := x1*w2 - w1*x2 + y1*z2 - z1*y2;

    # Omega_1 = x0^2 + x1^2 + x2^2 + x3^2
    rels[6]  := 1/3*(-x1*x2 + y1*y2 + z1*z2 + w1*w2) + 2/3*y1*y2 + 2/3*((1+alpha)/(1-beta))*z1*z2 + 2/3*((1-alpha)/(1+gamm))*w1*w2; ;


    rels[7]  := x2*y3 - y2*x3 - alpha*z2*w3 + alpha*w2*z3;
    rels[8]  := x2*z3 - z2*x3 - beta*w2*y3 + beta*y2*w3;
    rels[9] := x2*z3 - z2*x3 + w2*y3 - y2*w3;
    rels[10] := x2*w3 - w2*x3 + gamm*y2*z3 - gamm*z2*y3;
    rels[11] := x2*w3 - w2*x3 + y2*z3 - z2*y3;
    rels[12] := 1/3*(-x2*x3 + y2*y3 + z2*z3 + w2*w3) +2/3*y2*y3 + 2/3*((1+alpha)/(1-beta))*z2*z3 + 2/3*((1-alpha)/(1+gamm))*w2*w3;;

    I := TwoSidedIdeal(kQ, rels);
    A := kQ/I;
fi;
