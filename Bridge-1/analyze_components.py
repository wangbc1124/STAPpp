import re, csv, math
from collections import defaultdict

with open('Bridge-1.out', 'r') as f:
    content = f.read()
lines = content.split('\n')
stap_disps = {}
in_disp = False
for line in lines:
    if 'D I S P L A C E M E N T S' in line:
        in_disp = True; continue
    if in_disp:
        m = re.match(r'\s+(\d+)\s+(-?[\d.]+(?:e[+-]?\d+)?)\s+(-?[\d.]+(?:e[+-]?\d+)?)\s+(-?[\d.]+(?:e[+-]?\d+)?)\s+(-?[\d.]+(?:e[+-]?\d+)?)\s+(-?[\d.]+(?:e[+-]?\d+)?)\s+(-?[\d.]+(?:e[+-]?\d+)?)', line)
        if m:
            vals = [float(m.group(g)) for g in [2,3,4,5,6,7]]
            stap_disps[int(m.group(1))] = vals
        elif 'S T R E S S' in line or 'T O T A L' in line:
            break

stap_coords = {}
with open('Bridge-1.dat', 'r') as f:
    for line in f:
        p = line.strip().split()
        if len(p) == 10 and p[0].isdigit():
            nid = int(p[0])
            stap_coords[nid] = (float(p[7]), float(p[8]), float(p[9]))

abaqus_coords = {}
with open('../abaqus/all_instances_nodes.csv', 'r') as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
        abaqus_coords[(row[0], int(row[1]))] = (float(row[2]), float(row[3]), float(row[4]))

abaqus_disps = {}
with open('../abaqus/odb_displacements.csv', 'r') as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
        abaqus_disps[(row[1], int(row[0]))] = (float(row[2]), float(row[3]), float(row[4]))

stap_rev = defaultdict(list)
for nid, (x,y,z) in stap_coords.items():
    stap_rev[(round(x,4), round(y,4), round(z,4))].append(nid)

abaqus_to_stap = {}
for (inst, nid), (ax, ay, az) in abaqus_coords.items():
    key = (round(ax,4), round(ay,4), round(az,4))
    if key in stap_rev:
        abaqus_to_stap[(inst, nid)] = stap_rev[key][0]

# Per-component breakdown
comp_map = {}
for (inst, nid) in abaqus_to_stap:
    if 'FLOOR' in inst:
        comp_map[abaqus_to_stap[(inst, nid)]] = 'Deck-S4R'
    elif 'PIER' in inst:
        comp_map[abaqus_to_stap[(inst, nid)]] = 'Pier-C3D8R'
    elif 'SUPPORTBEAM' in inst:
        comp_map[abaqus_to_stap[(inst, nid)]] = 'Beam-B31'
    elif 'CABLE' in inst:
        comp_map[abaqus_to_stap[(inst, nid)]] = 'Cable-T3D2'
    elif 'RIVERBANK' in inst:
        comp_map[abaqus_to_stap[(inst, nid)]] = 'RiverBank-C3D8R'

comp_errs = defaultdict(lambda: {'n':0, 'a_dz':[], 's_dz':[], 'a_dy':[], 's_dy':[], 'a_dx':[], 's_dx':[]})

for (inst, nid), snid in abaqus_to_stap.items():
    if (inst, nid) not in abaqus_disps or snid not in stap_disps:
        continue
    if snid not in comp_map:
        continue
    comp = comp_map[snid]
    a_dx, a_dy, a_dz = abaqus_disps[(inst, nid)]
    s_dx, s_dy, s_dz, _, _, _ = stap_disps[snid]
    d = comp_errs[comp]
    d['n'] += 1
    d['a_dz'].append(a_dz); d['s_dz'].append(s_dz)
    d['a_dy'].append(a_dy); d['s_dy'].append(s_dy)
    d['a_dx'].append(a_dx); d['s_dx'].append(s_dx)

print('=' * 100)
print('  各部件位移误差分析')
print('=' * 100)
for comp in sorted(comp_errs.keys()):
    d = comp_errs[comp]
    a_dz_abs = [abs(v) for v in d['a_dz']]
    s_dz_abs = [abs(v) for v in d['s_dz']]
    a_dy_abs = [abs(v) for v in d['a_dy']]
    s_dy_abs = [abs(v) for v in d['s_dy']]

    dz_ratio = sum(s_dz_abs)/(sum(a_dz_abs)+1e-20)
    dy_ratio = sum(s_dy_abs)/(sum(a_dy_abs)+1e-20)

    # R^2 for dz
    mean_a = sum(d['a_dz'])/d['n']
    mean_s = sum(d['s_dz'])/d['n']
    ss_tot = sum((x-mean_s)**2 for x in d['s_dz'])
    ss_res = sum((d['s_dz'][i]-d['a_dz'][i])**2 for i in range(d['n']))
    r2_dz = 1 - ss_res/(ss_tot+1e-20)

    # Linear regression slope (S = k * A)
    num = sum(d['a_dz'][i]*d['s_dz'][i] for i in range(d['n']))
    den = sum(x*x for x in d['a_dz']) + 1e-20
    slope_dz = num/den

    # Same for dy
    num_y = sum(d['a_dy'][i]*d['s_dy'][i] for i in range(d['n']))
    den_y = sum(x*x for x in d['a_dy']) + 1e-20
    slope_dy = num_y/den_y

    print(f'\n{comp} (n={d["n"]}):')
    print(f'  DZ: sum|S|/sum|A|={dz_ratio:.3f}  slope={slope_dz:.3f}  R2={r2_dz:.3f}  '
          f'max|A|={max(a_dz_abs):.4e}  max|S|={max(s_dz_abs):.4e}')
    print(f'  DY: sum|S|/sum|A|={dy_ratio:.3f}  slope={slope_dy:.3f}  '
          f'max|A|={max(a_dy_abs):.4e}  max|S|={max(s_dy_abs):.4e}')

    # Per-node relative error distribution
    rel_errs = []
    for i in range(d['n']):
        a_mag = math.sqrt(d['a_dx'][i]**2 + d['a_dy'][i]**2 + d['a_dz'][i]**2)
        s_mag = math.sqrt(d['s_dx'][i]**2 + d['s_dy'][i]**2 + d['s_dz'][i]**2)
        if a_mag > 1e-6:
            rel_errs.append(abs(s_mag-a_mag)/a_mag)
    rel_errs.sort()
    if rel_errs:
        print(f'  相对误差分布: median={rel_errs[len(rel_errs)//2]:.3f}  '
              f'p25={rel_errs[len(rel_errs)//4]:.3f}  '
              f'p75={rel_errs[len(rel_errs)*3//4]:.3f}  '
              f'p90={rel_errs[int(len(rel_errs)*0.9)]:.3f}')

# Stress comparison
print('\n' + '=' * 100)
print('  应力对比')
print('=' * 100)
parts = re.split(r'S T R E S S  C A L C U L A T I O N S  F O R  E L E M E N T  G R O U P\s+(\d+)', content)
stap_stresses = {}
for i in range(1, len(parts), 2):
    grp = int(parts[i])
    text = parts[i+1]
    if grp in (1, 2):
        for line in text.split('\n'):
            m = re.match(r'\s+(\d+)\s+([-\d.e+]+)\s+([-\d.e+]+)', line)
            if m:
                stap_stresses[(grp, int(m.group(1)))] = {'type':'bar','force':float(m.group(2)),'stress':float(m.group(3))}
    elif grp in (4, 5):
        lines2 = text.split('\n'); j = 2
        while j < len(lines2):
            m1 = re.match(r'\s+(\d+)\s+([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)', lines2[j])
            if m1 and j+1 < len(lines2):
                m2 = re.match(r'\s+([-\d.e+]+)\s+([-\d.e+]+)\s+([-\d.e+]+)', lines2[j+1])
                if m2:
                    sx=float(m1.group(2));sy=float(m1.group(3));sz=float(m1.group(4))
                    txy=float(m2.group(1));tyz=float(m2.group(2));tzx=float(m2.group(3))
                    vm=math.sqrt(0.5*((sx-sy)**2+(sy-sz)**2+(sz-sx)**2+6*(txy**2+tyz**2+tzx**2)))
                    stap_stresses[(grp, int(m1.group(1)))]={'type':'h8','sx':sx,'sy':sy,'sz':sz,'txy':txy,'tyz':tyz,'tzx':tzx,'vm':vm}
                    j+=2; continue
            j+=1

abaqus_stresses = []
with open('../abaqus/odb_stresses.csv', 'r') as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
        abaqus_stresses.append((row[2], int(row[1]), float(row[9])))

abaq_inst_vm = defaultdict(lambda: {'max':0, 'sum':0, 'n':0})
for inst, eid, vm in abaqus_stresses:
    aba = abaq_inst_vm[inst]
    aba['n'] += 1; aba['sum'] += vm
    if vm > aba['max']: aba['max'] = vm

inst_to_grp = {'FLOOR':3, 'PIER':(4,5), 'RIVERBANK':5, 'CABLE':2, 'SUPPORTBEAM':1}
for inst in sorted(abaq_inst_vm.keys()):
    aba = abaq_inst_vm[inst]
    grp = None
    for key in inst_to_grp:
        if key in inst:
            grp = inst_to_grp[key]
            break
    if grp:
        if isinstance(grp, tuple):
            s_vms = [v['vm'] for k,v in stap_stresses.items() if k[0] in grp]
        else:
            s_vms = [v['vm'] for k,v in stap_stresses.items() if k[0] == grp]
        s_max = max(s_vms)/1e6 if s_vms else 0
        a_max = aba['max']/1e6
        ratio = s_max/(a_max+1e-20)
        print(f'{inst:<35s} Abaqus={a_max:.2f} MPa  STAP++={s_max:.2f} MPa  ratio={ratio:.3f}')

print('\nDone.')
