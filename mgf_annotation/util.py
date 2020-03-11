

from pyteomics.mass import nist_mass


def calculate_Mass_Diff():
    c_diff = nist_mass['C'][13][0] - nist_mass['C'][12][0]
    n_diff = nist_mass['N'][15][0] - nist_mass['N'][14][0]
    return(4 * c_diff + n_diff)


def calculate_allowed_Mass_Diff(n=1, z=1):
    m_diff = calculate_Mass_Diff()
    out = {}

    for i in range(1, n + 1):
        out[i] = {}
        for j in range(1, z + 1):
                out[i][j] = m_diff * i / j

    return(out)


def mass_diff_decharging_stuff(n=1, z=1):
    t = calculate_allowed_Mass_Diff(n=n, z=z)
    d = {0: {'combos': [{'z': 1, 'n': -1}], 'unique': True, 'decharge': {'state': False, 'z': 1}}}
    for n in t:
        for z in t[n]:
            if t[n][z] not in d:
                d[t[n][z]] = {'combos': [{'z': z, 'n': n}]}
            else:
                d[t[n][z]]['combos'].append({'z': z, 'n': n})

    for i in d:
        d[i]['unique'] = True if len(d[i]['combos']) == 1 else False
        if d[i]['unique']:
            d[i]['decharge'] = {'state': True, 'z': d[i]['combos'][0]['z']} if d[i]['combos'][0]['z'] > 1 else {'state': False, 'z': -1}
        else:
            d[i]['decharge'] = {'state': False, 'z': -1}
    return d


def gen_allowed_mass_diff(n=1, z=1):
    t = calculate_allowed_Mass_Diff(n=n, z=z)
    for i in t.keys():
        for j in t[i].keys():
            yield t[i][j]


def gen_allowed_mass_diff_with_sign(n=1, z=1):
    for i in gen_allowed_mass_diff(n, z):
        for sign in ['plus', 'minus']:
            if sign == 'plus':
                yield i
            else:
                yield -1 * i


def parse_scan_id(string):
    if 'Scan' in string:
        return string.split(' ')[2]
    else:
        return string.replace("\"", "").split('=')[-1]
