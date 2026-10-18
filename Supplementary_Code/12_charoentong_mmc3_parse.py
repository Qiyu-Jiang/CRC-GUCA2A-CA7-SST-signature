# -*- coding: utf-8 -*-
"""Parse mmc3 Table S6 -> Charoentong 28 immune cell gene sets -> CSV/JSON."""
import openpyxl, json, csv
from collections import OrderedDict

p = r".\outputs/phase1\09.charoentong\mmc3.xlsx"
wb = openpyxl.load_workbook(p, read_only=True)
ws = wb["Sheet1"]
rows = list(ws.iter_rows(values_only=True))
sets = OrderedDict()
immunity = {}
for r0 in rows[3:]:
    if r0[0] is None or r0[1] is None:
        continue
    gene = str(r0[0]).strip().upper()
    ct = str(r1 if False else r0[1]).strip()
    imm = str(r0[2]).strip() if r0[2] else ""
    sets.setdefault(ct, [])
    if gene not in sets[ct]:
        sets[ct].append(gene)
    immunity[ct] = imm
wb.close()

print("cell types:", len(sets))
for ct, genes in sets.items():
    print(f"  {ct:32s} {len(genes):3d} genes  [{immunity.get(ct,'')}]")

with open(r".\outputs/phase1\02.ssgsea\charoentong28_sets.json", "w", encoding="utf-8") as f:
    json.dump({ct: {"genes": gs, "immunity": immunity.get(ct, "")} for ct, gs in sets.items()},
              f, ensure_ascii=False, indent=1)
with open(r".\outputs/phase1\02.ssgsea\charoentong28.gmt", "w", encoding="utf-8") as f:
    for ct, gs in sets.items():
        f.write(ct + "\t" + immunity.get(ct, "") + "\t" + "\t".join(gs) + "\n")
print("saved: charoentong28_sets.json / charoentong28.gmt")
