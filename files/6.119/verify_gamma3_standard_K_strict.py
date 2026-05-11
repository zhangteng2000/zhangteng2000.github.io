import csv
from collections import defaultdict
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path
import mpmath as mp

iv = mp.iv
iv.dps = 80
mp.mp.dps = 80
getcontext().prec = 120

BASE_DIR = Path(__file__).resolve().parent
G_CERT = BASE_DIR / 'certificate_G_intervals_strict.csv'
FINAL_CERT = BASE_DIR / 'certificate_final_strip_strict.csv'

R_LEFT = Decimal('0.50414')
R_MID = Decimal('0.5')
R_RIGHT = Decimal('1')

# Use the exact algebraic value H(r0) = 16/81 * (28 + 19*sqrt(19))
# as an interval quantity, rather than a floating-point approximation.
H0 = iv.mpf(16) / 81 * (28 + 19 * iv.sqrt(19))

def g_dec(r: Decimal) -> Decimal:
    return Decimal(1) - r * r + Decimal(2) * (Decimal(1) - r) ** 3

def g_iv(r):
    return 1 - r * r + 2 * (1 - r) ** 3

def E0(r, b):
    k = iv.sqrt(4 - 3 * b * b)
    return (
        12 - iv.mpf(32) / 3 * r * r + iv.mpf(4) / 3 * b * b + iv.mpf(8) / 3 * b * r
    ) * (iv.mpf(8) / 3 + iv.mpf(3) / 2 * k - iv.mpf(35) / 12 * r) - 64 * b * r

def U1(r, b):
    k = iv.sqrt(4 - 3 * b * b)
    return (
        12 * (1 - r * r)
        + iv.mpf(4) / 3 * (r * r + 2 * b * r + b * b)
        + 4 * (1 - b * b) * k
        + k / 4 * (r * r + 6 * b * r + 9 * b * b)
        + iv.mpf(8) / 3 * r * (r * r + iv.mpf(7) / 16 * b * b - iv.mpf(35) / 32 * b * r)
        + iv.mpf(23) / 4 * b**3
    )

def Z(r, b):
    k = iv.sqrt(4 - 3 * b * b)
    return b * r * (iv.mpf(8) / 3 + iv.mpf(3) / 2 * k - iv.mpf(35) / 12 * r)

def alpha(r, b):
    return 12 - iv.mpf(32) / 3 * r * r + iv.mpf(4) / 3 * b * b

def beta(r, b):
    return iv.mpf(8) / 3 * b * r

def Hrbt(r, b, t):
    Zrb = Z(r, b)
    bet = beta(r, b)
    return -alpha(r, b) * t + iv.sqrt((Zrb - bet * t) ** 2 + 64 * b * b * r * r * t * (2 - t))

def F(r, b, u):
    t = (1 - r) * u
    return U1(r, b) - Z(r, b) + Hrbt(r, b, t) - H0

def verify_G():
    with G_CERT.open(newline='') as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise SystemExit('empty G certificate')

    cur = R_LEFT
    min_lower = None
    for row in rows:
        a = Decimal(row['r_left'])
        b = Decimal(row['r_right'])
        if a != cur:
            raise SystemExit(f'G coverage failed: expected interval to start at {cur}, got {a}')
        cur = b

        r = iv.mpf([str(a), str(b)])
        val = E0(r, g_iv(r))
        if val.a <= 0:
            raise SystemExit(f'G-certificate failed on {row}')
        if min_lower is None or val.a < min_lower:
            min_lower = val.a
    if cur != R_RIGHT:
        raise SystemExit(f'G coverage failed: final right endpoint is {cur}, not {R_RIGHT}')
    return min_lower

def verify_gap_patch():
    # Supplementary exact check for the analytically treated short band
    # 0.5 <= r <= 0.50414:
    # E0(r,1) = (280 r^3 - 470 r^2 - 826 r + 500)/9
    # is strictly decreasing on [0,1], so it suffices to check r = 0.50414.
    #
    # Use the exact rational identity 0.50414 = 25207 / 50000.
    r = Fraction(25207, 50000)
    val = (280 * r**3 - 470 * r**2 - 826 * r + 500) / 9
    if val <= 0:
        raise SystemExit('gap patch failed at r = 0.50414 for E0(r,1)')
    return val

def midpoint(a: Decimal, b: Decimal) -> Decimal:
    return (a + b) / 2

def verify_group_coverage(r_left: Decimal, r_right: Decimal, rows):
    b0 = g_dec(r_right)
    b_points = {b0, Decimal('1')}
    u_points = {Decimal('0'), Decimal('1')}
    boxes = []
    for row in rows:
        bl = Decimal(row['b_left'])
        br = Decimal(row['b_right'])
        ul = Decimal(row['u_left'])
        ur = Decimal(row['u_right'])
        if bl < b0 or br > Decimal('1') or ul < Decimal('0') or ur > Decimal('1'):
            raise SystemExit(f'box escapes its band rectangle: {row}')
        boxes.append((bl, br, ul, ur))
        b_points.add(bl)
        b_points.add(br)
        u_points.add(ul)
        u_points.add(ur)

    b_list = sorted(b_points)
    u_list = sorted(u_points)

    for i in range(len(b_list) - 1):
        left = max(b_list[i], b0)
        right = b_list[i + 1]
        if right <= left:
            continue
        bm = midpoint(left, right)
        for j in range(len(u_list) - 1):
            down = max(u_list[j], Decimal('0'))
            up = min(u_list[j + 1], Decimal('1'))
            if up <= down:
                continue
            um = midpoint(down, up)
            covered = False
            for bl, br, ul, ur in boxes:
                if bl <= bm <= br and ul <= um <= ur:
                    covered = True
                    break
            if not covered:
                raise SystemExit(
                    f'coverage failed in r-band [{r_left},{r_right}] '
                    f'at cell b in [{left},{right}], u in [{down},{up}]'
                )

def verify_final_strip():
    with FINAL_CERT.open(newline='') as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise SystemExit('empty final-strip certificate')

    groups = defaultdict(list)
    max_upper = None
    for row in rows:
        key = (Decimal(row['r_left']), Decimal(row['r_right']))
        groups[key].append(row)

        r = iv.mpf([row['r_left'], row['r_right']])
        b = iv.mpf([row['b_left'], row['b_right']])
        u = iv.mpf([row['u_left'], row['u_right']])
        val = F(r, b, u)
        if val.b >= 0:
            raise SystemExit(f'final-strip certificate failed on {row}')
        if max_upper is None or val.b > max_upper:
            max_upper = val.b

    # Check that the r-bands are contiguous and cover [0.50414, 1].
    band_keys = sorted(groups)
    cur = R_LEFT
    for a, b in band_keys:
        if a != cur:
            raise SystemExit(f'r-band coverage failed: expected start {cur}, got {a}')
        cur = b
    if cur != R_RIGHT:
        raise SystemExit(f'r-band coverage failed: final right endpoint is {cur}, not {R_RIGHT}')

    # Check, for each r-band, that the union of listed b-u rectangles covers [g(r_right), 1] x [0, 1].
    for key in band_keys:
        verify_group_coverage(key[0], key[1], groups[key])

    return max_upper, len(rows), len(groups)

if __name__ == '__main__':
    min_lower_G = verify_G()
    min_lower_gap = verify_gap_patch()
    max_upper_F, nboxes, nbands = verify_final_strip()
    print('All interval certificates verified.')
    print(f'G-interval minimum lower endpoint: {min_lower_G}')
    print(f'Gap-patch exact rational value for E0(0.50414,1): {min_lower_gap}')
    print(f'Final-strip maximum upper endpoint: {max_upper_F}')
    print(f'Final-strip boxes: {nboxes}; r-bands: {nbands}')
