import sympy as sp
from Output_cleaning import read_gap_file,fix_powers_in_file,replace_exp_with_zeta,fix_gap_file
from fractions import Fraction
from AllRels import rels,typ
# from A_ZALGEBRAS import rels, typ

req_basis, full_basis, mod_spec = read_gap_file(f"Codes/{typ}/{typ}")
# print(rels)

basis_len2 = [p for p in full_basis if len(p) == 4]
basis_len3 = [p for p in full_basis if len(p) == 6]
basis_len2_req = [p for p in req_basis if len(p) == 4]
basis_len3_req = [p for p in req_basis if len(p) == 6]

vertices=full_basis[:4]
arrows=full_basis[4:4*4]

one_hot_vertices = {
    v: [1 if i == j else 0 for j in range(len(vertices))]
    for i, v in enumerate(vertices)
}

arrows_1to2=arrows[:4]
arrows_2to3=arrows[4:8]
arrows_3to4=arrows[8:]


arrows_1to2_enc = {
    v: [1 if i == j else 0 for j in range(len(vertices))]
    for i, v in enumerate(arrows_1to2)
}
arrows_2to3_enc = {
    v: [1 if i == j else 0 for j in range(len(vertices))]
    for i, v in enumerate(arrows_2to3)
}
arrows_3to4_enc = {
    v: [1 if i == j else 0 for j in range(len(vertices))]
    for i, v in enumerate(arrows_3to4)
}


paths_1to3 = [p for p in basis_len2 if p[1] == '1' and p[3] == '2']
paths_2to4 = [p for p in basis_len2 if p[1] == '2' and p[3] == '3']


paths_1to3_req = [p for p in basis_len2_req if p[1] == '1' and p[3] == '2']
paths_2to4_req = [p for p in basis_len2_req if p[1] == '2' and p[3] == '3']

paths_1to2to3to4 = [p for p in basis_len3 if p[1] == '1' and p[3] == '2' and p[5] == '3']
paths_1to2to3to4_req = [p for p in basis_len3_req if p[1] == '1' and p[3] == '2' and p[5] == '3']



path_index = {p: i for i, p in enumerate(basis_len2)}

def rel_vector(rel):
    v = sp.zeros(1, len(basis_len2))
    for c, p in rel:
        v[0, path_index[p]] += c
    return v


R = sp.Matrix.vstack(*[rel_vector(r) for r in rels])

R_rref, pivots = R.rref()
# print(R_rref)
def path_type(path):
    # examples: x1y2 → 1→3, w2x3 → 2→4
    if path[1] == "1":
        return 1
    elif path[1] == "2":
        return 2
    else:
        raise ValueError(f"Unknown path type: {path}")

def encode_path_exact(path):
    # one-hot in full len2 space
    v = sp.zeros(len(basis_len2), 1)
    v[path_index[path], 0] = 1

    # build B
    B = sp.zeros(len(basis_len2), len(basis_len2_req))
    for j, b in enumerate(basis_len2_req):
        B[path_index[b], j] = 1

    # solve v = B c + R^T α
    M = B.row_join(R.T)
    sol = M.gauss_jordan_solve(v)[0]

    # --- NEW PART: slice by path type ---
    pt = path_type(path)
    type1 = len(paths_1to3_req)

    if pt == 1:
        return sp.Matrix(sol[:type1])
    elif pt == 2:
        return sp.Matrix(sol[type1:type1 + len(paths_2to4_req)]).T
    else:
        raise ValueError(f"Unknown path type for {path}")

new_rels=[]
for rel in rels:
    for a in arrows:
      if int(a[1])<int(rel[0][1][1]):
        new_rel=[]
        for r in rel:
          new_rel.append((r[0],a+r[1]))
        new_rels.append(new_rel)
      if int(a[1])>int(rel[0][1][3]):
        new_rel=[]
        for r in rel:
          new_rel.append((r[0],r[1]+a))
        new_rels.append(new_rel)

path_index3 = {p: i for i, p in enumerate(basis_len3)}

def rel_vector(rel):
    v = sp.zeros(1, len(basis_len3))
    for c, p in rel:
        v[0, path_index3[p]] += c
    return v


R3 = sp.Matrix.vstack(*[rel_vector(r) for r in new_rels])
# print(R3)

R_rref3, pivots3 = R3.rref()
print("here")

def encode_path_exact_3(path):
    v = sp.zeros(len(basis_len3), 1)
    v[path_index3[path], 0] = 1

    # Solve v = B c + R^T α exactly
    B = sp.zeros(len(basis_len3), len(basis_len3_req))
    for j, b in enumerate(basis_len3_req):
        B[path_index3[b], j] = 1

    M = B.row_join(R3.T)
    sol = M.gauss_jordan_solve(v)[0]

    return sp.Matrix(sol[:len(basis_len3_req)]).T

mod_spec_final={}
for mod in mod_spec:
  mod_spec_final[mod]=[]
  for m in mod_spec[mod]:
    # print("M:",len(m))
    if m in vertices:
      mod_spec_final[mod].append(list(one_hot_vertices[m]))
    elif m in arrows_1to2:
      mod_spec_final[mod].append(list(arrows_1to2_enc[m]))
    elif m in arrows_2to3:
      mod_spec_final[mod].append(list(arrows_2to3_enc[m]))
    elif m in arrows_3to4:
      mod_spec_final[mod].append(list(arrows_3to4_enc[m]))
    elif m in paths_1to3 or  m in paths_2to4:
      mod_spec_final[mod].append(list(encode_path_exact(m)))
    elif m in paths_1to2to3to4:
      mod_spec_final[mod].append(list(encode_path_exact_3(m)))
# print(mod_spec_final)

final_coef=[]
for mod in mod_spec_final:
  l=[mod,mod_spec_final[mod]]
  final_coef.append(l)

def write_gap(obj):
    """Recursively convert Python object to GAP code string"""
    if isinstance(obj, list):
        return "[ " + ", ".join(write_gap(x) for x in obj) + " ]"
    elif isinstance(obj, dict):
        items = ", ".join(f"{write_gap(k)} => {write_gap(v)}" for k, v in obj.items())
        return f"[ {items} ]"
    elif isinstance(obj, str):
        return f'"{obj}"'  # wrap Python string in GAP quotes
    else:
        return str(obj)  # numbers remain as-is

with open(f"Codes/{typ}/final_coef.g", "w") as f:
    f.write("final_coef := " + write_gap(final_coef) + ";\n")
replace_exp_with_zeta(f"Codes/{typ}/final_coef.g", n=3)
fix_powers_in_file(f"Codes/{typ}/final_coef.g")
fix_gap_file(f"Codes/{typ}/final_coef.g")
print("Written successfully.")

