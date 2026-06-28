import sympy as sp
import re

def normalize(expr):
    """Convert '(1)*x1*x2*x3' -> 'x1x2x3'"""
    if not isinstance(expr, str):
        raise TypeError(f"Expected string, got {type(expr)}: {expr}")
    return expr.replace("(1)*", "").replace("*", "").strip()

# def dump_raw(filename, basis, full_basis, modspec):
#     with open(filename, "w") as f:
#         f.write("basis = ")
#         f.write(repr(basis))
#         f.write("\n\n")

#         f.write("full_basis = ")
#         f.write(repr(full_basis))
#         f.write("\n\n")

#         f.write("modspec = ")
#         f.write(repr(modspec))
#         f.write("\n")

#     print(f"Raw data written to {filename}")




import re
def read_gap_file(filename):
    with open(filename, "r") as f:
        text = f.read()

    # ---------------- BASIS ----------------
    basis_match = re.search(r"Basis:\[(.*?)\]\s*MOD_SPEC:", text, re.S)
    if not basis_match:
        raise ValueError("Could not find Basis block")

    basis_raw = basis_match.group(1)
    basis_elems = re.findall(r"\[\(1\)\*([^\]]+)\]", basis_raw)
    basis = [normalize(b) for b in basis_elems]

    # ---------------- FULL BASIS ----------------
    full_basis_match = re.search(r"Full Basis:\[(.*?)\]", text, re.S)
    if full_basis_match:
        full_basis_raw = full_basis_match.group(1)
        full_basis_elems = re.findall(r"\(\d+\)\*([a-zA-Z0-9_]+(?:\*[a-zA-Z0-9_]+)*)", full_basis_raw)
        full_basis = [normalize(b) for b in full_basis_elems]
    else:
        full_basis = []

    # ---------------- MOD_SPEC ----------------
    modspec = {}
    mod_match = re.search(r"MOD_SPEC:\[(.*)\]\s*$", text, re.S)
    if not mod_match:
        raise ValueError("Could not find MOD_SPEC block")

    mod_raw = mod_match.group(1)
    entries = re.findall(r'\[\s*"([^"]+)"\s*,\s*\[(.*?)\]\s*\]', mod_raw, re.S)
    for key, value_block in entries:
        elems = re.findall(r"\(\d+\)\*([a-zA-Z0-9_]+(?:\*[a-zA-Z0-9_]+)*)", value_block)
        modspec[key] = [normalize(e) for e in elems]

    return basis, full_basis, modspec


# ========== USE ==========
# basis, full_basis, modspec = read_gap_file("4VerticesCommutativeRels/4VerticesCommutativeRels")
# dump_raw("4VerticesCommutativeRels/parsed_output.txt", basis, full_basis, modspec)

# print("Basis size:", len(basis))
# print("Full Basis size:", len(full_basis))
# print("Sample Full Basis elements:", full_basis[:10])
# print(basis)
# print("\n\n",full_basis,"\n\n",modspec)

def fix_powers_in_file(filepath):
    with open(filepath, "r") as f:
        content = f.read()
    
    if "**" in content:
        content = content.replace("**", "^")
        with open(filepath, "w") as f:
            f.write(content)
        print("Fixed ** -> ^ in", filepath)
    else:
        print("No ** found in", filepath)

# Usage

def replace_exp_with_zeta(filepath, n=9):
    """Cleanup pass — catch any remaining exp() or E(n) forms and fix double negatives"""
    with open(filepath, "r") as f:
        content = f.read()
    # Replace any remaining E(n)^k
    content = re.sub(rf"E\({n}\)\^(\d+)", r"E(\1)", content)  # temporary
    content = re.sub(rf"E\({n}\)", f"E({n})", content)         # restore

    # Fix double negatives
    content = content.replace("--", "")

    # Replace remaining exp() forms
    replacements = {}
    for k in range(n):
        zeta_k = sp.exp(2 * sp.pi * sp.I * k / n)
        gap_str = f"E({n})^{k}" if k > 1 else (f"E({n})" if k == 1 else "1")

        for candidate in [str(zeta_k), str(sp.simplify(zeta_k)), str(sp.nsimplify(zeta_k))]:
            replacements[candidate] = gap_str
        for candidate in [str(-zeta_k), str(sp.simplify(-zeta_k))]:
            replacements[candidate] = f"-{gap_str}"

    for old, new in sorted(replacements.items(), key=lambda x: -len(x[0])):
        if old in content:
            content = content.replace(old, new)

    with open(filepath, "w") as f:
        f.write(content)
    print(f"Done: {filepath}")
        

def fix_gap_file(filepath):
    with open(filepath, "r") as f:
        content = f.read()
    
    # Fix double negatives
    content = content.replace("--", "")
    
    # Fix powers
    content = content.replace("**", "^")
    
    # Remove trailing backslashes and extra whitespace lines
    lines = content.splitlines()
    cleaned_lines = [line for line in lines if line.strip() not in ["\\", ""]]
    content = "\n".join(cleaned_lines) + "\n"

    with open(filepath, "w") as f:
        f.write(content)
    print("Fixed file:", filepath)
    # print(f"Done replacing exp -> zeta in {filepath}")
