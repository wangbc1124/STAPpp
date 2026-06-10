import sys
with open('Bridge-1-3d-fix.dat','r') as f:
    lines = f.readlines()
nodes = {}
for i in range(2, 4087):
    p = lines[i].strip().split()
    if len(p) >= 10:
        try:
            nodes[int(p[0])] = (float(p[7]), float(p[8]), float(p[9]))
        except:
            pass
bad = 0
for j in range(4090, min(4110, len(lines))):
    p = lines[j].strip().split()
    if len(p) >= 4:
        try:
            n1 = int(p[1])
            n2 = int(p[2])
            p1 = nodes.get(n1, (0, 0, 0))
            p2 = nodes.get(n2, (0, 0, 0))
            dy = abs(p1[1] - p2[1])
            if dy > 0.01:
                bad += 1
                print("BAD Elem %s: N%d(%.0f,%.4f,%.0f) -> N%d(%.0f,%.4f,%.0f) dY=%.4f" % (
                    p[0], n1, p1[0], p1[1], p1[2], n2, p2[0], p2[1], p2[2], dy))
        except:
            pass
if bad == 0:
    print("All 20 cables OK")
else:
    print("%d cables BAD" % bad)
