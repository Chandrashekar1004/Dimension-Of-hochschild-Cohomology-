LoadPackage("qpa");
Read("Codes/Functions.g");

start_time := Runtime();

Q:=MakeQuiver(4);

# Read("Codes/A-ZALGEBRAS.g");
Read("Codes/ALGEBRAS.g");


main_path:="/gap-4.14.0/Codes/";

mod_spec := ComputeModuleSpec( A, type, main_path );

Exec("/myenv/bin/python3 /gap-4.14.0/Codes/Coefficients.py");

Read(Concatenation("/gap-4.14.0/Codes/",type,"/","final_coef.g"));

final_coef:=CoeffConversion(final_coef,K);;

proj_res:=ComputeProjRes(A,K,final_coef);;

ans1:=ComputeCohomology(proj_res,1);
ans2:=ComputeCohomology(proj_res,2);
ans3:=ComputeCohomology(proj_res,3);
# ans4:=ComputeCohomology(proj_res,4);


end_time := Runtime();
elapsed := end_time - start_time;

output := OutputTextFile(Concatenation("gap-4.14.0/Codes/", type, "/cohomology.txt"), false);;
AppendTo(output, "dim1 = ", ans1, "\n");
AppendTo(output, "dim2 = ", ans2, "\n");
AppendTo(output, "dim3 = ", ans3, "\n");


AppendTo(output, "runtime = ", elapsed, " ms\n");
CloseStream(output);
