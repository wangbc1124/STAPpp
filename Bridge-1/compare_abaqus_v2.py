"""Abaqus ODB vs STAP++-3D comprehensive comparison."""
import csv, re, math, sys
from collections import defaultdict

OUT_PATH = r'D:\STAPpp-master\STAPpp-master\Bridge-1\Bridge-1-3d-fix.out'
DAT_PATH = r'D:\STAPpp-master\STAPpp-master\Bridge-1\Bridge-1-3d-fix.dat'
ABAQUS_DIR = r'D:\STAPpp-master\STAPpp-master\abaqus'

# ========== 1. Load STAP++ displacements ==========
print("=== Loading STAP++ data ===")
with open(OUT_PATH, 'r') as f:
    content = f.read()

stap_disps = {}
lines = content.split('\n')
in_disp = False
for line in lines:
    if 'D I S P L A C E M E N T S' in line:
        in_disp = True; continue
    if in_disp:
        m = re.match(r'\s+(\d+)\s+(-?[\d.]+(?:e[+-]?\d+)?)\s+(-?[\d.]+(?:e[+-]?\d+)?)\s+(-?[\d.]+(?:e[+-]?\d+)?)\s+(-?[\d.]+(?:e[+-]?\d+)?)\s+(-?[\d.]+(?:e[+-]?\d+)?)\s+(-?[\d.]+(?:e[+-]?\d+)?)', line)
        if m:
            vals = [float(m.group(g)) for g in [2,3,4,5,6,7]]
            stap_disps[int(m.group(1))] = tuple(vals)
        elif 'S T R E S S' in line or 'T O T A L' in line:
            break

print(f"  STAP++ displacements: {len(stap_disps)} nodes")

# STAP++ stresses by element group
parts = re.split(r'S T R E S S  C A L C U L A T I O N S  F O R  E L E M E N T  G R O U P\s+(\d+)', content)
stap_stresses = {}
for i in range(1, len(parts), 2):
    grp = int(parts[i])
    text = parts[i+1]
    if grp in (1, 2):  # Bar
        for line in text.split('\n'):
            m = re.match(r'\s+(\d+)\s+([-\d.e+]+)\s+([-\d.e+]+)', line)
            if m:
                eid = int(m.group(1)); f_val = float(m.group(2)); s_val = float(m.group(3))
                stap_stresses[(grp, eid)] = {'type':'bar','force':f_val,'stress':s_val,'vm':abs(s_val)}
    elif grp == 3:  # Q4
        for line in text.split('\n'):
            m = re.match(r'\s+(\d+)\s+([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)', line)
            if m:
                sx=float(m.group(2)); sy=float(m.group(3)); txy=float(m.group(4))
                vm=math.sqrt(sx**2+sy**2-sx*sy+3*txy**2)
                stap_stresses[(grp, int(m.group(1)))]={'type':'q4','sx':sx,'sy':sy,'txy':txy,'vm':vm}
    elif grp in (4, 5):  # H8
        lines2 = text.split('\n'); j = 2
        while j < len(lines2):
            m1 = re.match(r'\s+(\d+)\s+([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)', lines2[j])
            if m1 and j+1 < len(lines2):
                m2 = re.match(r'\s+([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)', lines2[j+1])
                if m2:
                    sx=float(m1.group(2)); sy=float(m1.group(3)); sz=float(m1.group(4))
                    txy=float(m2.group(1)); tyz=float(m2.group(2)); tzx=float(m2.group(3))
                    vm=math.sqrt(0.5*((sx-sy)**2+(sy-sz)**2+(sz-sx)**2+6*(txy**2+tyz**2+tzx**2)))
                    stap_stresses[(grp, int(m1.group(1)))]={'type':'h8','sx':sx,'sy':sy,'sz':sz,'txy':txy,'tyz':tyz,'tzx':tzx,'vm':vm}
                    j+=2; continue
            j+=1
    elif grp in (9,):  # Beam3D
        for line in text.split('\n'):
            m = re.match(r'\s+(\d+)\s+([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)', line)
            if m and 'SHELL' not in text and 'ELEMENT' not in text and 'SIGMA' not in text:
                pass
            # Beam3D stress format: 8 values on one line

print(f"  STAP++ stress entries: {len(stap_stresses)}")

# STAP++ node coordinates
stap_coords = {}
with open(DAT_PATH, 'r') as f:
    dat_lines = f.readlines()
for line in dat_lines[2:]:
    p = line.strip().split()
    if len(p) == 10:
        nid = int(p[0])
        stap_coords[nid] = (float(p[7]), float(p[8]), float(p[9]))

print(f"  STAP++ coordinates: {len(stap_coords)} nodes")

# ========== 2. Load Abaqus ODB data ==========
print("\n=== Loading Abaqus data ===")

# Node coordinates from CSV (maps to global instances)
abaqus_coords = {}
with open(f'{ABAQUS_DIR}/all_instances_nodes.csv', 'r') as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
        inst = row[0]; nid = int(row[1])
        abaqus_coords[(inst, nid)] = (float(row[2]), float(row[3]), float(row[4]))

# Abaqus displacements
abaqus_disps = {}
with open(f'{ABAQUS_DIR}/odb_displacements.csv', 'r') as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
        abaqus_disps[(row[1], int(row[0]))] = (float(row[2]), float(row[3]), float(row[4]))

# Abaqus stresses (element integration points, partial)
abaqus_stresses = []
with open(f'{ABAQUS_DIR}/odb_stresses.csv', 'r') as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
        abaqus_stresses.append((row[2], int(row[1]), float(row[3]), float(row[4]),
                                 float(row[5]), float(row[6]), float(row[7]),
                                 float(row[8]), float(row[9])))

print(f"  Abaqus coords: {len(abaqus_coords)}, disp: {len(abaqus_disps)}, stress: {len(abaqus_stresses)}")

# ========== 3. Map Abaqus -> STAP++ nodes ==========
print("\n=== Mapping nodes ===")
# Build STAP++ reverse coordinate index
stap_rev = defaultdict(list)
for nid, (x,y,z) in stap_coords.items():
    stap_rev[(round(x,4), round(y,4), round(z,4))].append(nid)

abaqus_to_stap = {}
unmapped = 0
for (inst, nid), (ax, ay, az) in abaqus_coords.items():
    key = (round(ax,4), round(ay,4), round(az,4))
    if key in stap_rev:
        abaqus_to_stap[(inst, nid)] = stap_rev[key][0]
    else:
        unmapped += 1

print(f"  Mapped: {len(abaqus_to_stap)}, Unmapped: {unmapped}")

# ========== 4. Displacement comparison by instance ==========
print("\n" + "="*90)
print("  位移对比: Abaqus ODB vs STAP++-3D")
print("="*90)

per_inst = defaultdict(lambda: {'n': 0, 'a_dz_max': 0, 's_dz_max': 0, 'diffs_dz': [], 'diffs_dy': [], 'diffs_dx': []})
all_diffs = []

for (inst, nid), snid in abaqus_to_stap.items():
    if (inst, nid) not in abaqus_disps or snid not in stap_disps:
        continue
    a_dx, a_dy, a_dz = abaqus_disps[(inst, nid)]
    s_dx, s_dy, s_dz, *_ = stap_disps[snid]
    d = per_inst[inst]
    d['n'] += 1
    adz, sdz = abs(a_dz), abs(s_dz)
    if adz > d['a_dz_max']: d['a_dz_max'] = adz
    if sdz > d['s_dz_max']: d['s_dz_max'] = sdz
    d['diffs_dz'].append(abs(a_dz - s_dz))
    d['diffs_dy'].append(abs(a_dy - s_dy))
    d['diffs_dx'].append(abs(a_dx - s_dx))
    all_diffs.append((inst, nid, snid, a_dz, s_dz, a_dx, s_dx, a_dy, s_dy))

# Summary table
inst_names = {
    'PART-FLOOR-1': '甲板', 'PART-PIER-1': '桥塔1', 'PART-PIER-2': '桥塔2',
    'PART-SUPPORTBEAM-1': '支撑梁1', 'PART-SUPPORTBEAM-2': '支撑梁2',
    'PART-RIVERBANK-1': '河岸1', 'PART-RIVERBANK-2': '河岸2',
}
cable_insts = [k for k in per_inst if 'CABLE' in k]

print(f"\n{'部件':<25s} {'节点数':>6s} {'Abaqus Max|DZ|':>15s} {'STAP++ Max|DZ|':>15s} {'比值':>8s} {'Avg|diff_dz|':>13s}")
print("-" * 85)
for inst in sorted(per_inst.keys()):
    d = per_inst[inst]
    ratio = d['s_dz_max']/d['a_dz_max'] if d['a_dz_max'] > 1e-15 else 0
    avg_diff = sum(d['diffs_dz'])/len(d['diffs_dz']) if d['diffs_dz'] else 0
    label = inst_names.get(inst, inst)
    print(f"{label:<25s} {d['n']:>6d} {d['a_dz_max']:>15.6e} {d['s_dz_max']:>15.6e} {ratio:>8.4f} {avg_diff:>13.6e}")

# Cable summary
cable_data = [(inst, per_inst[inst]) for inst in cable_insts if inst in per_inst]
if cable_data:
    cable_n = sum(d['n'] for _, d in cable_data)
    cable_a_max = max(d['a_dz_max'] for _, d in cable_data)
    cable_s_max = max(d['s_dz_max'] for _, d in cable_data)
    ratio = cable_s_max/cable_a_max if cable_a_max > 1e-15 else 0
    print(f"{'索(Cable) x20':<25s} {cable_n:>6d} {cable_a_max:>15.6e} {cable_s_max:>15.6e} {ratio:>8.4f}")

# Global stats
all_adz = [abs(a_dz) for _,_,_,a_dz,_,_,_,_,_ in all_diffs]
all_sdz = [abs(s_dz) for _,_,_,_,s_dz,_,_,_,_ in all_diffs]
all_adz_sum = sum(all_adz)
all_sdz_sum = sum(all_sdz)
print(f"\n  全局总 |DZ|: Abaqus={all_adz_sum:.4e}, STAP++={all_sdz_sum:.4e}, 比值={all_sdz_sum/all_adz_sum if all_adz_sum>0 else 0:.4f}")
print(f"  全局 Max |DZ|: Abaqus={max(all_adz):.4e}, STAP++={max(all_sdz):.4e}")

# ========== 5. Node-level errors ==========
print("\n" + "="*90)
print("  节点级误差分析 (Top 20 |dz_abaqus - dz_stap|)")
print("="*90)

all_diffs.sort(key=lambda x: abs(x[3]-x[4]), reverse=True)
print(f"\n{'Abaqus(inst,node)':<32s} {'STAP++':>6s} {'A_dz(m)':>12s} {'S_dz(m)':>12s} {'Diff(m)':>12s}")
print("-" * 80)
for inst, nid, snid, adz, sdz, *_ in all_diffs[:20]:
    print(f"  {inst}:{str(nid):<20s} {snid:>6d} {adz:>12.6e} {sdz:>12.6e} {adz-sdz:>12.6e}")

# ========== 6. Key node comparison ==========
print("\n" + "="*90)
print("  关键节点对比")
print("="*90)

# Find key nodes by searching for specific Abaqus coords
key_targets = {
    '甲板中心': (0, 10, 0),
    '甲板中心(-10)': (0, -10, 0),
    '甲板(-225,10)': (-225, 10, 0),
    '甲板(225,10)': (225, 10, 0),
    '塔1顶(Z=150)': (0, -10, 150),
    '塔2顶(Z=150)': (0, 10, 150),
    '塔1底(Z=-50)': (0, -10, -50),
}

for desc, target in key_targets.items():
    tx, ty, tz = target
    best_abaq = None
    best_dist = 1e9
    for (inst, nid), (ax, ay, az) in abaqus_coords.items():
        d = math.sqrt((ax-tx)**2 + (ay-ty)**2 + (az-tz)**2)
        if d < best_dist:
            best_dist = d
            best_abaq = (inst, nid, ax, ay, az)
    if best_abaq and best_dist < 1.0:
        inst, nid, ax, ay, az = best_abaq
        if (inst, nid) in abaqus_disps and (inst, nid) in abaqus_to_stap:
            snid = abaqus_to_stap[(inst, nid)]
            a_dz = abaqus_disps[(inst, nid)][2]
            s_dz = stap_disps.get(snid, (0,0,0,0,0,0))[2]
            a_dy = abaqus_disps[(inst, nid)][1]
            s_dy = stap_disps.get(snid, (0,0,0,0,0,0))[1]
            print(f"{desc}: Abaqus dz={a_dz:.4e}, STAP++ dz={s_dz:.4e}, 比值={abs(s_dz/(a_dz+1e-30)):.4f}")
            print(f"         Abaqus dy={a_dy:.4e}, STAP++ dy={s_dy:.4e}")

# ========== 7. Stress comparison ==========
print("\n" + "="*90)
print("  应力对比")
print("="*90)

# Per-instance Abaqus stress
abaq_inst_vm = defaultdict(lambda: {'max':0, 'sum':0, 'n':0})
for inst, eid, s11, s22, s33, s12, s23, s13, vm in abaqus_stresses:
    aba = abaq_inst_vm[inst]
    aba['n'] += 1; aba['sum'] += vm
    if vm > aba['max']: aba['max'] = vm

print(f"\n{'部件':<25s} {'Abaqus数':>8s} {'Abaqus max VM(MPa)':>19s} {'STAP++ max VM(MPa)':>19s} {'类型':>10s}")
print("-" * 85)
for inst in sorted(abaq_inst_vm.keys()):
    aba = abaq_inst_vm[inst]
    label = inst_names.get(inst, inst)
    # Map to STAP++ group
    if 'FLOOR' in inst:
        s_max = max(v['vm'] for k,v in stap_stresses.items() if k[0]==3)/1e6 if any(k[0]==3 for k in stap_stresses) else 0
        stype = 'Q4'
    elif 'PIER' in inst:
        s_max = max(v['vm'] for k,v in stap_stresses.items() if k[0] in (4,5))/1e6 if any(k[0] in (4,5) for k in stap_stresses) else 0
        stype = 'H8'
    elif 'RIVERBANK' in inst:
        s_max = max(v['vm'] for k,v in stap_stresses.items() if k[0]==5)/1e6 if any(k[0]==5 for k in stap_stresses) else 0
        stype = 'H8'
    elif 'SUPPORTBEAM' in inst:
        s_max = max(v['vm'] for k,v in stap_stresses.items() if k[0]==1)/1e6 if any(k[0]==1 for k in stap_stresses) else 0
        stype = 'Bar'
    elif 'CABLE' in inst:
        s_max = max(v['vm'] for k,v in stap_stresses.items() if k[0]==2)/1e6 if any(k[0]==2 for k in stap_stresses) else 0
        stype = 'Bar'
    else:
        s_max = 0; stype = '?'
    print(f"{label:<25s} {aba['n']:>8d} {aba['max']/1e6:>19.4f} {s_max:>19.4f} {stype:>10s}")

# Print overall von Mises stats
print(f"\n  Abaqus 全局 max VM: {max(d['max'] for d in aba_inst_vm.values())/1e6:.2f} MPa")
stap_vms = [v['vm'] for v in stap_stresses.values()]
print(f"  STAP++ 全局 max VM: {max(stap_vms)/1e6:.2f} MPa")

print("\n=== Done ===")
