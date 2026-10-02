#!/usr/bin/env python3
"""
build_ortholog_map.py — MGI HOM_MouseHumanSequence.rpt → 小鼠→人类基因符号映射表

输入: MGI 同源表（DB Class Key 聚合同源组；每行一个物种成员）
输出: mouse_to_human_symbol.tsv
  列: mouse_symbol, human_symbol, is_1to1
  - 同源组内恰好 1 鼠 + 1 人 → is_1to1=True（主分析使用）
  - 多对多/多对一 → is_1to1=False（备用，默认不用于 scDRS 映射）
"""
import sys
from collections import defaultdict

SRC = sys.argv[1] if len(sys.argv) > 1 else "project/datasets/ortholog/HOM_MouseHumanSequence.rpt"
OUT = sys.argv[2] if len(sys.argv) > 2 else "project/datasets/ortholog/mouse_to_human_symbol.tsv"

groups = defaultdict(lambda: {"mouse": set(), "human": set()})
with open(SRC, encoding="utf-8", errors="replace") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    i_key = header.index("DB Class Key")
    i_org = header.index("Common Organism Name")
    i_sym = header.index("Symbol")
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) <= max(i_key, i_org, i_sym):
            continue
        key, org, sym = parts[i_key], parts[i_org], parts[i_sym]
        if not sym:
            continue
        if org.startswith("mouse"):
            groups[key]["mouse"].add(sym)
        elif org == "human":
            groups[key]["human"].add(sym)

n_pairs = n_1to1 = 0
seen_mouse = set()
with open(OUT, "w", encoding="utf-8") as out:
    out.write("mouse_symbol\thuman_symbol\tis_1to1\n")
    for key, g in groups.items():
        if not g["mouse"] or not g["human"]:
            continue
        one2one = len(g["mouse"]) == 1 and len(g["human"]) == 1
        for m in sorted(g["mouse"]):
            for h in sorted(g["human"]):
                out.write(f"{m}\t{h}\t{one2one}\n")
                n_pairs += 1
                if one2one:
                    n_1to1 += 1
                seen_mouse.add(m)

print(f"homology groups: {len(groups)}")
print(f"pairs written: {n_pairs} (1:1 pairs: {n_1to1})")
print(f"unique mouse symbols covered: {len(seen_mouse)}")
print(f"output: {OUT}")
