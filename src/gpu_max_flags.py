#!/usr/bin/env python3
"""Emit GPU_MAX_* preprocessor flags derived from the active input deck.

The Fortran code reads UCNS3D.DAT at runtime, but fixed GPU work-array
sizes must be selected at compile time. This helper mirrors the small part
of parameters.f90 that derives dimensions used by those arrays. With no
argument it checks the current build directory for UCNS3D.DAT and companion
files such as REALGAS.DAT.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

NUM_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[dDeE][-+]?\d+)?")

DEFAULTS = {
    "dimensiona": 3,
    "governingequations": 1,
    "turbulence": 1,
    "passivescalar": 0,
    "spatiladiscret": 3,
    "spatialorder": 5,
    "ees": 0,
    "nof_species": 0,
    "mp_modelc": 0,
}


def as_int(token: str) -> int:
    return int(float(token.replace("D", "E").replace("d", "e")))


def numeric_records(path: Path) -> list[list[str]]:
    records: list[list[str]] = []
    if not path.exists():
        return records
    for line in path.read_text(errors="ignore").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if not re.match(r"^[+-]?(?:\d|\.)", stripped):
            continue
        tokens = NUM_RE.findall(stripped)
        if tokens:
            records.append(tokens)
    return records


def first_numeric_after_skips(path: Path, skips: int) -> list[str] | None:
    if not path.exists():
        return None
    lines = path.read_text(errors="ignore").splitlines()[skips:]
    for line in lines:
        stripped = line.strip()
        if re.match(r"^[+-]?(?:\d|\.)", stripped):
            tokens = NUM_RE.findall(stripped)
            if tokens:
                return tokens
    return None


def derive(config: Path) -> dict[str, int]:
    values = dict(DEFAULTS)
    records = numeric_records(config)

    if len(records) > 0 and len(records[0]) >= 1:
        values["dimensiona"] = as_int(records[0][0])
    if len(records) > 1 and len(records[1]) >= 1:
        values["governingequations"] = as_int(records[1][0])
    if len(records) > 2:
        if len(records[2]) >= 1:
            values["turbulence"] = as_int(records[2][0])
        if len(records[2]) >= 3:
            values["passivescalar"] = as_int(records[2][2])
    if len(records) > 6:
        if len(records[6]) >= 1:
            values["spatiladiscret"] = as_int(records[6][0])
        if len(records[6]) >= 3:
            values["spatialorder"] = as_int(records[6][2])
    if len(records) > 7 and len(records[7]) >= 2:
        values["ees"] = as_int(records[7][1])

    base_dir = config.parent if config.parent != Path("") else Path(".")
    multispecies = (base_dir / "MULTISPECIES.DAT").exists()
    realgas = (base_dir / "REALGAS.DAT").exists()

    if multispecies:
        species_tokens = first_numeric_after_skips(base_dir / "MULTISPECIES.DAT", 2)
        if species_tokens:
            values["nof_species"] = as_int(species_tokens[0])
        model_tokens = first_numeric_after_skips(base_dir / "MULTISPECIES_DIFF.DAT", 1)
        if model_tokens:
            values["mp_modelc"] = as_int(model_tokens[0])

    if realgas:
        species_tokens = first_numeric_after_skips(base_dir / "REALGAS.DAT", 2)
        if species_tokens:
            values["nof_species"] = as_int(species_tokens[0])

    dim = values["dimensiona"]
    # 2D cases still compile shared 3D helper routines, so the fixed-size
    # non-xpu work arrays must keep 3D-capable storage even when dimensiona=2.
    storage_dim = max(3, dim)
    species = values["nof_species"]
    if values["governingequations"] <= 2:
        nof_variables = 5 if dim == 3 else 4
    else:
        nof_variables = 1

    if multispecies:
        flow_vars = 5 if dim == 3 else 4
        if values["mp_modelc"] == 0:
            nof_variables = flow_vars + species + max(species - 1, 0)
        else:
            nof_variables = flow_vars + species

    if realgas:
        nof_variables = dim + 3 + species

    iorder = max(1, values["spatialorder"] - 1)
    if values["spatiladiscret"] == 1 and values["spatialorder"] == 1:
        iorder = 1

    typesten = 7
    if values["spatiladiscret"] == 3:
        ees = values["ees"]
        if storage_dim == 2:
            if ees in (1, 2):
                typesten = 9
            else:
                typesten = 5
        else:
            if ees in (1, 2):
                typesten = 15
            else:
                typesten = 7

    # UCNS3D.DAT only toggles turbulence on/off; the concrete model is set in
    # parameters.f90 profiles. Reserve room for both SA (1) and k-omega SST (2).
    turbulence_equations = 2 if values["turbulence"] == 1 else 0

    # parameters.f90 sets igqrules from iorder, then gaussianpoints/quadalloc
    # translate that into numberofpoints and numberofpoints2 maxima.
    dg = 1
    igqrules = min(iorder + 1, 6) if dg == 1 else min(iorder, 6)
    qp3 = {
        1: (1, 1, 1, 1, 1, 1),
        2: (8, 4, 5, 6, 4, 3),
        3: (27, 14, 15, 18, 9, 6),
        4: (64, 24, 15, 24, 16, 10),
        5: (125, 35, 15, 50, 25, 15),
        6: (216, 56, 15, 60, 36, 21),
    }
    qp2 = {
        1: (1, 1, 1),
        2: (4, 3, 2),
        3: (9, 6, 3),
        4: (16, 10, 4),
        5: (25, 15, 5),
        6: (36, 21, 6),
    }
    if storage_dim == 3:
        qp_hexa, qp_tetra, qp_pyra, qp_prism, qp_quad, qp_triangle = qp3.get(igqrules, (216, 56, 15, 60, 36, 36))
        numberofpoints = max(qp_hexa, qp_tetra * 6, qp_pyra, qp_prism) if dg == 1 else max(qp_hexa, qp_tetra, qp_pyra, qp_prism)
        numberofpoints2 = max(qp_quad, qp_triangle)
        qp_all = max(qp_hexa, qp_tetra, qp_pyra, qp_prism, qp_quad, qp_triangle)
    else:
        qp_quad, qp_triangle, qp_line = qp2.get(igqrules, (36, 36, 9))
        numberofpoints = max(qp_quad, qp_triangle * 2) if dg == 1 else max(qp_quad, qp_triangle)
        numberofpoints2 = qp_line
        qp_all = max(qp_quad, qp_triangle, qp_line)

    extf = 3
    if storage_dim == 3:
        dg_dof = ((iorder + 1) * (iorder + 2) * (iorder + 3)) // 6
        idegfree = dg_dof - 1
        neighbours = 9 if iorder == 1 else dg_dof * extf
        neighbours2 = 7 if iorder == 1 else (12 if iorder <= 3 else 20)
    else:
        dg_dof = ((iorder + 1) * (iorder + 2)) // 2
        idegfree = dg_dof - 1
        neighbours = 5 if iorder == 1 else dg_dof * extf
        neighbours2 = 6 if iorder == 1 else (7 if iorder <= 3 else 11)

    return {
        "GPU_MAX_DIM": storage_dim,
        "GPU_MAX_IORDER": iorder,
        "GPU_MAX_TYPESTEN": typesten,
        "GPU_MAX_NVAR": max(5, nof_variables),
        "GPU_MAX_SPECIES": max(1, species),
        "GPU_MAX_TURBULENCE": max(1, turbulence_equations),
        "GPU_MAX_PASSIVE": max(1, values["passivescalar"]),
        "GPU_MAX_EXTF": extf,
        "GPU_MAX_IDEGFREE": max(1, idegfree),
        "GPU_MAX_DOF": max(1, dg_dof),
        "GPU_MAX_NEIGHBOURS": max(1, neighbours),
        "GPU_MAX_QP_FACE": max(1, numberofpoints2),
        "GPU_MAX_QP_VOLUME": max(1, numberofpoints),
        "GPU_MAX_QP_ALL": max(1, qp_all),
        "GPU_MAX_NEIGHBOURS2": neighbours2,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "config",
        nargs="?",
        default="UCNS3D.DAT",
        help="UCNS3D.DAT path; defaults to ./UCNS3D.DAT in the current build directory",
    )
    parser.add_argument("--summary", action="store_true")
    args = parser.parse_args()

    config = Path(args.config)
    vals = derive(config)
    if args.summary:
        for key in sorted(vals):
            print(f"{key}={vals[key]}")
    else:
        sys.stdout.write(" ".join(f"-D{key}={value}" for key, value in vals.items()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
