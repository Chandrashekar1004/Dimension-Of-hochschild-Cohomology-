# Read("Codes_26_02/GAP/Quiver.g");
# Read("Codes_26_02/GAP/AlgebraA/Parameters.g");
LoadPackage("qpa");

MakeQuiver := function( number_of_vertices )
    local vertices, arrows, Q;

    if number_of_vertices = 2 then
        Q := Quiver( ["v1", "v2"], [
            ["v1", "v2", "x1"],
            ["v1", "v2", "x2"]
        ] );

    elif number_of_vertices = 3 then
        Q := Quiver( ["v1", "v2", "v3"], [
            ["v1", "v2", "x1"],
            ["v1", "v2", "y1"],
            ["v1", "v2", "z1"],
            ["v2", "v3", "x2"],
            ["v2", "v3", "y2"],
            ["v2", "v3", "z2"]
        ] );

    elif number_of_vertices = 4 then
        Q := Quiver( ["v1", "v2", "v3", "v4"], [
            ["v1", "v2", "w1"],
            ["v1", "v2", "x1"],
            ["v1", "v2", "y1"],
            ["v1", "v2", "z1"],
            ["v2", "v3", "w2"],
            ["v2", "v3", "x2"],
            ["v2", "v3", "y2"],
            ["v2", "v3", "z2"],
            ["v3", "v4", "w3"],
            ["v3", "v4", "x3"],
            ["v3", "v4", "y3"],
            ["v3", "v4", "z3"]
        ] );

    else
        Error( "number_of_vertices must be 2, 3, or 4" );
    fi;

    return Q;
end;


ComputeModuleSpec := function( A, type, main_path )

    local output, env, Q1, QxQ, A1, basis1, basis_vectors_1, basis, basis_vectors,
          vertices, vertex_indices, vector_spaces, i, j, dimension_vector,
          arrows, module_specification, a, coeffs, Q,dir_path,output_path,
          make_map;


    dir_path := Concatenation(main_path, type);

    if not IsDirectoryPath( dir_path ) then
        CreateDir( dir_path );
    fi;

    output_path:=Concatenation(dir_path,"/",type);
    output := OutputTextFile( output_path, true );;
    env := EnvelopingAlgebra( A );
    Q   := QuiverOfPathAlgebra( A );
    Q1  := QuiverOfPathAlgebra( A );
    QxQ := QuiverOfPathAlgebra( env );
    A1  := PathAlgebra( Rationals, Q1 );

    basis1          := CanonicalBasis( A1 );
    basis_vectors_1 := BasisVectors( basis1 );
    basis           := CanonicalBasis( A );
    basis_vectors   := BasisVectors( basis );

    AppendTo( output, "Basis:", basis_vectors, "\n" );
    AppendTo( output, "Full Basis:", basis_vectors_1, "\n" );

    vertices       := VerticesOfQuiver( Q );
    vertex_indices := [ 1 .. Length( vertices ) ];
    vector_spaces  := NullMat( Length( vertices ), Length( vertices ) );

    for i in vertex_indices do
        for j in vertex_indices do
            vector_spaces[ i ][ j ] := PositionsNonzero( vertices[ i ] * basis_vectors * vertices[ j ] );
        od;
    od;

    dimension_vector := List( Concatenation( vector_spaces ), Length );

    make_map := function( a )
        local components, source, target, source_i, target_i,
              source_space, target_space, dims, map_on_basis;

        components := [ ProjectFromProductQuiver( 1, a ),
                        ProjectFromProductQuiver( 2, a ) ];
        source   := List( components, SourceOfPath );
        target   := List( components, TargetOfPath );
        source_i := List( source, VertexIndex );
        target_i := List( target, VertexIndex );

        source_space := vector_spaces[ source_i[ 1 ] ][ source_i[ 2 ] ];
        target_space := vector_spaces[ target_i[ 1 ] ][ target_i[ 2 ] ];

        dims := [ Length( source_space ), Length( target_space ) ];
        if dims[ 1 ] = 0 or dims[ 2 ] = 0 then
            return fail;
        fi;

        map_on_basis := OppositePath( components[ 1 ] ) * basis_vectors{ source_space } * components[ 2 ];
        return map_on_basis;
    end;

    arrows             := ArrowsOfQuiver( QxQ );
    module_specification := [ ];

    for a in arrows do
        coeffs := [];
        if make_map( a ) <> fail then
            Add( module_specification, [ String( a ), make_map( a ) ] );
        fi;
    od;

    AppendTo( output, "\nMOD_SPEC:", module_specification );
    CloseStream( output );

    return module_specification;

end;


CoeffConversion := function( final_coef, k )
    local one, i, j, m;
    
    one := One( k );
    
    for i in [1..Length(final_coef)] do
        for j in [1..Length(final_coef[i][2])] do
            for m in [1..Length(final_coef[i][2][j])] do
                final_coef[i][2][j][m] := one * final_coef[i][2][j][m];
            od;
        od;
    od;
    
    return final_coef;
end;

ComputeProjRes := function( A, k, final_coef )
    local env, basis, basis_vectors, vertices, vertex_indices,
          vector_spaces, i, j, dimension_vector,
          Amod, proj_res, Q;

    Q             := QuiverOfPathAlgebra( A );
    env           := EnvelopingAlgebra( A );
    basis         := CanonicalBasis( A );
    basis_vectors := BasisVectors( basis );
    vertices      := VerticesOfQuiver( Q );
    vertex_indices := [ 1 .. Length( vertices ) ];

    vector_spaces := NullMat( Length( vertices ), Length( vertices ) );

    for i in vertex_indices do
        for j in vertex_indices do
            vector_spaces[ i ][ j ] := PositionsNonzero( vertices[ i ] * basis_vectors * vertices[ j ] );
        od;
    od;

    dimension_vector := List( Concatenation( vector_spaces ), Length );

    final_coef := CoeffConversion( final_coef, k );

    Amod     := RightModuleOverPathAlgebra( env, dimension_vector, final_coef );
    proj_res := ProjectiveResolution( Amod );

    return proj_res;
end;

FlattenDeltaImages := function( delta_images )
    local del_mat, del, row, i, j, k;

    del_mat := [];

    for del in delta_images do
        row := [];
        for i in del do
            for j in i do
                for k in j do
                    Add( row, k );
                od;
            od;
        od;
        Add( del_mat, row );
    od;

    return del_mat;
end;


ComputeCohomology := function( proj_res, dim )
    local obj_A, P_n, P_nM1, P_nP1,
          H_n, H_nM1, H_nP1,
          d_n, d_n_mat, d_nP1, d_nP1_mat,
          delta_images, delta_images_nP1,
          del_mat, del_nP1_mat,
          H_n_mat, H_nP1_mat, H_nP1_flattened, H_n_flattened,
          LambanP1, Lamba,
          ker_dim, img_dim, ans,
          f, h, del;

    obj_A  := ObjectOfComplex( proj_res, -1 );
    P_n    := ObjectOfComplex( proj_res, dim );
    P_nM1  := ObjectOfComplex( proj_res, dim - 1 );
    P_nP1  := ObjectOfComplex( proj_res, dim + 1 );

    H_n   := HomOverAlgebra( P_n,   obj_A );
    H_nM1 := HomOverAlgebra( P_nM1, obj_A );
    H_nP1 := HomOverAlgebra( P_nP1, obj_A );

    d_n      := DifferentialOfComplex( proj_res, dim );
    d_n_mat  := MatricesOfPathAlgebraMatModuleHomomorphism( d_n );
    d_nP1     := DifferentialOfComplex( proj_res, dim + 1 );
    d_nP1_mat := MatricesOfPathAlgebraMatModuleHomomorphism( d_nP1 );

    # Compute delta images for H_nM1
    delta_images := [];
    for f in H_nM1 do
        del := d_n * f;
        Add( delta_images, MatricesOfPathAlgebraMatModuleHomomorphism( del ) );
    od;

    H_n_mat := [];
    for h in H_n do
        Add( H_n_mat, MatricesOfPathAlgebraMatModuleHomomorphism( h ) );
    od;

    # Compute kernel dimension
    if Length( H_nP1 ) <> 0 then
        delta_images_nP1 := [];
        for f in H_n do
            del := d_nP1 * f;
            Add( delta_images_nP1, MatricesOfPathAlgebraMatModuleHomomorphism( del ) );
        od;

        H_nP1_mat := [];
        for h in H_nP1 do
            Add( H_nP1_mat, MatricesOfPathAlgebraMatModuleHomomorphism( h ) );
        od;

        del_nP1_mat      := FlattenDeltaImages( delta_images_nP1 );
        H_nP1_flattened  := FlattenDeltaImages( H_nP1_mat );
        LambanP1         := List( del_nP1_mat, row -> SolutionMat( H_nP1_flattened, row ) );
        ker_dim          := Length( NullspaceMat( LambanP1 ) );
    else
        ker_dim := Length( H_n );
    fi;

    # Compute image dimension
    del_mat       := FlattenDeltaImages( delta_images );
    H_n_flattened := FlattenDeltaImages( H_n_mat );
    Lamba         := List( del_mat, row -> SolutionMat( H_n_flattened, row ) );
    img_dim       := RankMat( Lamba );

    ans := ker_dim - img_dim;

    return ans;
#     return rec(
#     ans          := ker_dim - img_dim,
#     H_n_mat      := H_n_mat,
#     H_n_flat     := H_n_flattened,
#     Lamba        := Lamba,
#     LambanP1     := LambanP1,
#     K            := NullspaceMat(LambanP1)
# );;
end;