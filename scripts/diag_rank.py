#!/usr/bin/env python3
# §862-диаг: сопоставление ранжиров схем нормировки с «истинным этажом»
# (измеренной ошибкой огрубления). Вход: build/diag/<сцена>_lp{0..3}.e (дампы
# per-piece E по лестнице lp) + <сцена>_diag.e (hits, sdelta за прогон).
import sys, re

def load(path):
    d = {}
    for ln in open(path):
        m = re.match(r"E1 p=(\d+) tri=(-?\d+) e=(\S+) hits=(\S+) sdelta=(\S+)", ln)
        if m:
            d[int(m.group(1))] = (float(m.group(3)), float(m.group(4)), float(m.group(5)))
    return d

def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            r[order[k]] = avg
        i = j + 1
    return r

def spearman(a, b):
    if len(a) < 3:
        return float("nan")
    ra, rb = ranks(a), ranks(b)
    n = len(a)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    da = sum((x - ma) ** 2 for x in ra) ** 0.5
    db = sum((y - mb) ** 2 for y in rb) ** 0.5
    return num / (da * db) if da > 0 and db > 0 else float("nan")

for scene in sys.argv[1:]:
    base = load(f"build/diag/{scene}_lp0.e")
    lad = {lp: load(f"build/diag/{scene}_lp{lp}.e") for lp in (1, 2, 3)}
    diag = load(f"build/diag/{scene}_diag.e")
    ps = sorted(base)
    maxE = max(base[p][0] for p in ps) if ps else 0.0
    sel, true_fl, vA, vB, vV, mats = [], [], [], [], [], []
    for p in ps:
        e0, hits, sd = base[p]
        if e0 <= 1e-9 * max(maxE, 1e-30):
            continue  # ненагруженный кусок: ему нечего терять, ранжир не определён
        fl = 0
        for lp in (1, 2, 3):
            if abs(lad[lp][p][0] - e0) <= 0.05 * e0:
                fl = lp
        d = diag[p]
        sel.append(p)
        true_fl.append(fl)
        vA.append(d[1])                          # (а) число событий
        vB.append(d[2] / d[1] if d[1] else 0.0)  # (б) среднее Δ за событие
        vV.append(e0)                            # (в) энергия (депозит)
        mats.append((p, e0, d[1], d[2]))
    print(f"== {scene}: кусков в ранжире {len(sel)} (из {len(ps)}), maxE={maxE:.4g}")
    tf = [0] * 4
    for f in true_fl:
        tf[f] += 1
    print(f"   истинные этажи: 0:{tf[0]} 1:{tf[1]} 2:{tf[2]} 3:{tf[3]}")
    for name, v in (("(а) события", vA), ("(б) среднее Δ", vB), ("(в) энергия", vV)):
        print(f"   {name:14s} Спирмен = {spearman(v, true_fl):+.3f}")
    # НК: перемешать
    import random
    random.seed(7)
    w = vA[:]
    random.shuffle(w)
    print(f"   {'НК перемеш.':14s} Спирмен = {spearman(w, true_fl):+.3f} (~0 ожидается)")
