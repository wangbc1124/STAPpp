import re, csv, math
from collections import defaultdict

# Load STAP++ nodes (exactly N_NODES)
with open('Bridge-1.dat', 'r') as f:
    lines = f.readlines()
N_NODES = int(lines[1].strip().split()[0])
stap_coords = {}
for i in range(2, 2 + N_NODES):
    p = lines[i].strip().split()
    if len(p) == 10 and p[0].isdigit():
        nid = int(p[0])
        stap_coords[nid] = (float(p[7]), float(p[8]), float(p[9]))

# Load STAP++ displacements
with open('Bridge-1.out', 'r') as f:
    content = f.read()
out_lines = content.split('\n')
stap_disps = {}
in_disp = False
for line in out_lines:
    if 'D I S P L A C E M E N T S' in line:
        in_disp = True; continue
    if in_disp:
        m = re.match(r'\s+(\d+)\s+(-?[\d.]+(?:e[+-]?\d+)?)\s+(-?[\d.]+(?:e[+-]?\d+)?)\s+(-?[\d.]+(?:e[+-]?\d+)?)\s+(-?[\d.]+(?:e[+-]?\d+)?)\s+(-?[\d.]+(?:e[+-]?\d+)?)\s+(-?[\d.]+(?:e[+-]?\d+)?)', line)
        if m:
            vals = [float(m.group(g)) for g in [2,3,4,5,6,7]]
            stap_disps[int(m.group(1))] = vals
        elif 'S T R E S S' in line or 'T O T A L' in line:
            break

# Load Abaqus data
abaqus_coords = {}
abaqus_disps = {}
with open('../abaqus/all_instances_nodes.csv', 'r') as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
        abaqus_coords[(row[0], int(row[1]))] = (float(row[2]), float(row[3]), float(row[4]))
with open('../abaqus/odb_displacements.csv', 'r') as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
        abaqus_disps[(row[1], int(row[0]))] = (float(row[2]), float(row[3]), float(row[4]))

# Map
stap_rev = defaultdict(list)
for nid, (x,y,z) in stap_coords.items():
    stap_rev[(round(x,4), round(y,4), round(z,4))].append(nid)

# FLOOR comparison
print('=' * 80)
print('  DECK (FLOOR / Shell4) DISPLACEMENT COMPARISON')
print('=' * 80)
floor_data = []
for (inst, nid), (ax, ay, az) in abaqus_coords.items():
    if 'FLOOR' not in inst:
        continue
    key = (round(ax,4), round(ay,4), round(az,4))
    if key not in stap_rev:
        continue
    snid = stap_rev[key][0]
    if (inst, nid) not in abaqus_disps or snid not in stap_disps:
        continue
    a_dx, a_dy, a_dz = abaqus_disps[(inst, nid)]
    s_dx, s_dy, s_dz, s_rx, s_ry, s_rz = stap_disps[snid]
    floor_data.append((ax, ay, a_dx, a_dy, a_dz, s_dx, s_dy, s_dz, s_rz, snid))

floor_data.sort(key=lambda x: abs(x[4]), reverse=True)
n = len(floor_data)

# Global stats
a_dzs = [abs(d[4]) for d in floor_data]
s_dzs = [abs(d[7]) for d in floor_data]
a_dys = [abs(d[3]) for d in floor_data]
s_dys = [abs(d[6]) for d in floor_data]
dz_ratio = sum(s_dzs)/(sum(a_dzs)+1e-20)
dy_ratio = sum(s_dys)/(sum(a_dys)+1e-20)
rms_dz = math.sqrt(sum((d[4]-d[7])**2 for d in floor_data)/n)
rms_dy = math.sqrt(sum((d[3]-d[6])**2 for d in floor_data)/n)

print(f'Nodes: {n}')
print(f'Abaqus max|DZ|: {max(a_dzs):.6e} m')
print(f'STAP++ max|DZ|: {max(s_dzs):.6e} m')
print(f'DZ ratio sum|S|/sum|A|: {dz_ratio:.4f}')
print(f'DY ratio sum|S|/sum|A|: {dy_ratio:.4f}')
print(f'RMS dz error: {rms_dz:.6e} m')
print(f'RMS dy error: {rms_dy:.6e} m')

# Per-node relative error
rel_errs = []
for d in floor_data:
    a_mag = math.sqrt(d[2]**2 + d[3]**2 + d[4]**2)
    s_mag = math.sqrt(d[5]**2 + d[6]**2 + d[7]**2)
    if a_mag > 1e-6:
        rel_errs.append(abs(s_mag - a_mag) / a_mag)
rel_errs.sort()
print(f'Displacement rel error: median={rel_errs[n//2]:.3f} p25={rel_errs[n//4]:.3f} p75={rel_errs[3*n//4]:.3f}')

# Top 15 |A_dz|
print(f'\nTop 15 Abaqus |dz|:')
print(f'{"X":>8s} {"Y":>8s} {"A_dz":>12s} {"S_dz":>12s} {"Err_dz":>12s} {"S/A":>8s}')
print('-' * 65)
for ax, ay, a_dx, a_dy, a_dz, s_dx, s_dy, s_dz, s_rz, snid in floor_data[:15]:
    ratio = abs(s_dz/(a_dz+1e-20))
    print(f'{ax:>8.1f} {ay:>8.1f} {a_dz:>12.6e} {s_dz:>12.6e} {abs(a_dz-s_dz):>12.6e} {ratio:>8.4f}')

# Centerline profile (Y=10)
centerline = [(d[0], d[4], d[7]) for d in floor_data if abs(d[1]-10) < 0.01]
centerline.sort(key=lambda x: x[0])
print(f'\nDZ along centerline Y=10:')
print(f'{"X":>8s} {"Abaqus":>12s} {"STAP++":>12s} {"Ratio":>8s}')
print('-' * 45)
for x, adz, sdz in centerline[::5]:
    ratio = abs(sdz/(adz+1e-20))
    print(f'{x:>8.1f} {adz:>12.6e} {sdz:>12.6e} {ratio:>8.4f}')

# Check: are deck Z DOFs fixed in STAP++?
print(f'\nDeck node BC check:')
for nid in [1, 2, 9, 10, 19, 20]:
    if nid in stap_coords:
        x, y, z = stap_coords[nid]
        # Find BC from .dat
        for i in range(2, 2+N_NODES):
            p = lines[i].strip().split()
            if len(p) == 10 and int(p[0]) == nid:
                bc = [int(p[j]) for j in range(1, 7)]
                disp = stap_disps.get(nid, (0,0,0,0,0,0))
                print(f'  Node {nid} ({x:.1f},{y:.1f},{z:.1f}): BC={bc} disp_dz={disp[2]:.6e}')
                break

# Check RZ constraint of deck nodes
rz_fixed = 0
rz_free = 0
for (inst, nid), (ax, ay, az) in abaqus_coords.items():
    if 'FLOOR' not in inst:
        continue
    key = (round(ax,4), round(ay,4), round(az,4))
    if key not in stap_rev:
        continue
    snid = stap_rev[key][0]
    for i in range(2, 2+N_NODES):
        p = lines[i].strip().split()
        if len(p) == 10 and int(p[0]) == snid:
            if int(p[6]) == 1:
                rz_fixed += 1
            else:
                rz_free += 1
            break
print(f'\nDeck RZ DOF: fixed={rz_fixed}, free={rz_free}')
