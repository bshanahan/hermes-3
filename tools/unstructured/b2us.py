"""Reader for GOAT / SOLPS unstructured grid files (traduit.out.b2us)."""
import numpy as np


def read_b2us(path):
    lines = open(path).read().split("\n")
    heads = [i for i, l in enumerate(lines) if l.startswith("*cf")]
    heads.append(len(lines))

    def body(i):
        j = heads[heads.index(i) + 1]
        return " ".join(lines[i + 1:j]).split()

    g = {}
    for i in heads[:-1]:
        h = lines[i]
        if "nCi,nFc" in h:
            g["nCi"], g["nFc"], g["nVx"], g["nCg"], g["nFs"], g["nFt"] = map(int, body(i))
        elif "isClassicalGrid" in h:
            g["classical"] = int(body(i)[0])
        elif h.startswith("*cf: Vx"):
            a = np.array(body(i), float).reshape(-1, 7)
            g["vx"], g["vy"], g["vpsi"], g["vdpsidx"], g["vdpsidy"], g["vffbz"] = a[:, 1:].T
        elif h.startswith("*cf: cv"):
            a = np.array(body(i), float)
            a = a.reshape(g["nCi"], -1)
            g["cvVxP"] = a[:, 1:3].astype(int)
            g["cx"], g["cy"], g["cpsi"], g["cbp"], g["cbt"] = a[:, 3:8].T
            g["cflags"] = a[:, 8:-2].astype(int)
            g["cReg"], g["cFt"] = a[:, -2].astype(int), a[:, -1].astype(int)
        elif h.startswith("*cf: fc"):
            a = np.array(body(i), int).reshape(-1, 6)
            g["fcVx"], g["fcLbl"], g["fcReg"], g["fcAligned"] = a[:, 1:3] - 1, a[:, 3], a[:, 4], a[:, 5]
        elif h.startswith("*cf: ft"):
            a = np.array(body(i), int).reshape(-1, 6)
            g["ftCvP"], g["ftFcP"], g["ftReg"] = a[:, 1:3], a[:, 3:5], a[:, 5]
        elif h.startswith("*cf: fs"):
            a = np.array(body(i), float).reshape(-1, 4)
            g["fsFcP"], g["fsPsi"] = a[:, 1:3].astype(int), a[:, 3]
        else:
            name = h.split()[-1]
            if name in ("cvVx", "cvFc", "ftCv", "ftFc", "fsFc"):
                g[name] = np.array(body(i), int) - 1
    start = g["cvVxP"][:, 0] - 1
    count = g["cvVxP"][:, 1]
    g["cells"] = [g["cvVx"][s:s + n] for s, n in zip(start, count)]
    g["cell_faces"] = [g["cvFc"][s:s + n] for s, n in zip(start, count)]
    g["tubes"] = [g["ftCv"][s - 1:s - 1 + n] for s, n in g["ftCvP"]]
    return g
