import csv
from collections import defaultdict

stap_coords = {}
with open('Bridge-1.dat', 'r') as f:
    lines = f.readlines()

# Find Shell4 group header
shell4_start = None
for i, line in enumerate(lines):
    if '10  400  1' in line:
        shell4_start = i
        break

# Read all nodes
for line in lines[2:]:
    p = line.strip().split()
    if len(p) == 10 and p[0].isdigit():
        nid = int(p[0])
        stap_coords[nid] = (float(p[7]), float(p[8]), float(p[9]))
    elif p and not p[0].isdigit():
        break

# Shell4 element connectivity
if shell4_start:
    print('Shell4 elements (first 10):')
    for i in range(shell4_start + 2, shell4_start + 12):
        p = lines[i].strip().split()
        if len(p) >= 5 and p[0].isdigit():
            eid = int(p[0])
            conn = [int(n) for n in p[1:5]]
            coords_str = ' '.join(f'{stap_coords.get(n, (0,0,0))}' for n in conn)
            print(f'  Elem {eid}: nodes={conn}')

shell4_nodes = set()
if shell4_start:
    for i in range(shell4_start + 2, shell4_start + 2 + 400):
        p = lines[i].strip().split()
        if len(p) >= 5 and p[0].isdigit():
            for n in p[1:5]:
                shell4_nodes.add(int(n))

print(f'\nUnique Shell4 nodes: {len(shell4_nodes)}')

# Coords of these nodes
xs = [stap_coords[n][0] for n in shell4_nodes]
ys = [stap_coords[n][1] for n in shell4_nodes]
zs = [stap_coords[n][2] for n in shell4_nodes]
print(f'X range: {min(xs):.1f} to {max(xs):.1f}')
print(f'Y range: {min(ys):.1f} to {max(ys):.1f}')
print(f'Z range: {min(zs):.4f} to {max(zs):.4f}')

# Sample Shell4 node coords
for nid in sorted(shell4_nodes)[:10]:
    c = stap_coords.get(nid, (0,0,0))
    print(f'  Node {nid}: ({c[0]:.2f}, {c[1]:.2f}, {c[2]:.2f})')

# THE CRITICAL CHECK: compare Abaqus FLOOR coords with Shell4 node coords
abaqus_coords = {}
with open('../abaqus/all_instances_nodes.csv', 'r') as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
        abaqus_coords[(row[0], int(row[1]))] = (float(row[2]), float(row[3]), float(row[4]))

shell4_coords = {n: stap_coords[n] for n in shell4_nodes}
print('\nMatching Abaqus FLOOR nodes to STAP++ Shell4 nodes:')
matched = 0
unmatched = 0
for (inst, nid), (ax, ay, az) in abaqus_coords.items():
    if 'FLOOR' not in inst:
        continue
    best_dist = 1e9
    best_snid = 0
    for snid, (sx, sy, sz) in shell4_coords.items():
        d = ((ax-sx)**2 + (ay-sy)**2 + (az-sz)**2)**0.5
        if d < best_dist:
            best_dist = d
            best_snid = snid
    if best_dist < 0.1:
        matched += 1
    else:
        unmatched += 1
        if unmatched <= 10:
            # Find closest ANY STAP++ node
            best_any = 1e9
            best_any_id = 0
            for snid, (sx, sy, sz) in stap_coords.items():
                d = ((ax-sx)**2 + (ay-sy)**2 + (az-sz)**2)**0.5
                if d < best_any:
                    best_any = d
                    best_any_id = snid
            sc = stap_coords.get(best_any_id, (0,0,0))
            print(f'  {inst}:{nid} ({ax:.2f},{ay:.2f},{az:.2f}) best_any={best_any_id} ({sc[0]:.2f},{sc[1]:.2f},{sc[2]:.2f}) dist={best_any:.4f}')

print(f'\nMatched: {matched}, Unmatched: {unmatched} out of {matched+unmatched} FLOOR nodes')
