from numpy import power, sqrt
from scipy.special import factorial


def wigner_3j(j1, j2, j, m1, m2, m):
    """Computes Wigner 3 j-symbols using
       https://mathworld.wolfram.com/Wigner3j-Symbol.html

    Args:
        j1: Top left 3 j-symbol value.
        j2: Top center 3 j-symbol value.
        j: Top right 3 j-symbol value.
        m1: Bottom left 3 j-symbol value.
        m2: Bottom center 3 j-symbol value.
        m: Bottom right 3 j-symbol value.

    Returns:
        Wigner 3 j-symbol value.
    """
    if not _selection_rules(j1, j2, j, m1, m2, m): return 0.
    try:
        if j1 == j2 and m1 == -1*m2 and j == 0 and m == 0:
            return _special_case1(j1, m1)
        elif m1 == 0 and m2 == 0 and m == 0:
            return _special_case2(j1, j2, j)
        elif j == j1 + j2:
            return _special_case3(j1, j2, m1, m2, m)
        elif j1 == j2 + j:
            if not _selection_rules(j2, j, j1, m2, m, m1): return 0.
            return _special_case3(j2, j, m2, m, m1)
        elif m1 == j1 and m2 == -j1:
            return _special_case4(j1, j2, j, m)
        elif m1 == j1 and m == -j1:
            if not _selection_rules(j1, j, j2, m1, m, m2): return 0.
            return power(-1., j1 + j2 + j)*_special_case4(j1, j, j2, m2)
        else:
            return _general_case(j1, j2, j, m1, m2, m)
    except OverflowError:
        print("{}  {}  {}\n{}  {}  {}".format(j1, j2, j, m1, m2, m))
        raise


def _selection_rules(j1, j2, j, m1, m2, m):
    """Checks selection rules.

    Args:
        j1: Top left 3 j-symbol value.
        j2: Top center 3 j-symbol value.
        j: Top right 3 j-symbol value.
        m1: Bottom left 3 j-symbol value.
        m2: Bottom center 3 j-symbol value.
        m: Bottom right 3 j-symbol value.

    Returns:
        True if values pass selection rules check, else False.
    """
    if not (m1 + m2 == -m or (abs(j1 - j2) <= j and j <= j1 + j2)):
        return False
    for x, y in zip([j1, j2, j], [m1, m2, -m]):
        if not (y >= -1*abs(x) and y <= abs(x)):
            return False
    return True


def _special_case1(l, m):
    """Calculates wigner 3-j symbol of the form:
       |l   l   0|
       |m  -m   0|

    Args:
        l: Top left 3 j-symbol value.
        m: Bottom left 3 j-symbol value.

    Returns:
        Wigner 3 j-symbol value.
    """
    return power(-1., l - m)/sqrt(2*l + 1)


def _special_case2(j1, j2, j):
    """Calculates wigner 3-j symbol of the form:
       |j1  j2  j|
       |0   0   0|

    Args:
        j1: Top left 3 j-symbol value.
        j2: Top center 3 j-symbol value.
        j: Top right 3 j-symbol value.

    Returns:
        Wigner 3 j-symbol value.
    """
    J = j1 + j2 + j
    if J % 2 == 0:
        g = J/2
        term1 = power(-1., g)
        term2 = sqrt((factorial(2*g - 2*j1, exact=True) *
                      factorial(2*g - 2*j2, exact=True) *
                      factorial(2*g - 2*j, exact=True)) /
                     factorial(2*g + 1, exact=True))
        term3 = factorial(g, exact=True) / \
            (factorial(g - j1, exact=True) *
             factorial(g - j2, exact=True) *
             factorial(g - j, exact=True))
        return term1*term2*term3
    else:
        return 0.


def _special_case3(j1, j2, m1, m2, m):
    """Calculates wigner 3-j symbol of the form:
       |j1  j2  j1+j2|
       |m1  m2    m  |

    Args:
        j1: Top left 3 j-symbol value.
        j2: Top center 3 j-symbol value.
        m1: Bottom left 3 j-symbol value.
        m2: Bottom center 3 j-symbol value.
        m: Bottom right 3 j-symbol value.

    Returns:
        Wigner 3 j-symbol value.
    """
    term1 = power(-1., j1 - j2 - m)
    term2 = sqrt((factorial(2*j1, exact=True)*factorial(2*j2, exact=True)) /
                 factorial(2*j1 + 2*j2 + 1, exact=True))
    term3 = sqrt((factorial(j1 + j2 - m, exact=True) *
                  factorial(j1 + j2 + m, exact=True)) /
                 (factorial(j1 + m1, exact=True) *
                  factorial(j1 - m1, exact=True) *
                  factorial(j2 + m2, exact=True) *
                  factorial(j2 - m2, exact=True)))
    return term1*term2*term3


def _special_case4(j1, j2, j, m):
    """Calculates wigner 3-j symbol of the form:
       |j1   j2  j|
       |j1  -j1  m|

    Args:
        j1: Top left 3 j-symbol value.
        j2: Top center 3 j-symbol value.
        j: Top right 3 j-symbol value.
        m: Bottom right 3 j-symbol value.

    Returns:
        Wigner 3 j-symbol value.
    """
    term1 = power(-1., -j1 + j2 - m)
    term2 = sqrt((factorial(2*j1, exact=True) *
                  factorial(-j1 + j2 + j, exact=True)) /
                 (factorial(j1 + j2 + j + 1, exact=True) *
                  factorial(j1 - j2 + j, exact=True)))
    term3 = sqrt((factorial(j1 + j2 - m, exact=True) *
                  factorial(j + m, exact=True)) /
                 (factorial(j1 + j2 - j, exact=True) *
                  factorial(-j1 + j2 + m, exact=True) *
                  factorial(j - m, exact=True)))
    return term1*term2*term3


def _general_case(a, b, c, alpha, beta, gamma):
    """Calculates wigner 3-j symbol for the general case:
       |  a     b      c  |
       |alpha  beta  gamma|

    Args:
        a: Top left 3 j-symbol value.
        b: Top center 3 j-symbol value.
        c: Top right 3 j-symbol value.
        alpha: Bottom left 3 j-symbol value.
        beta: Bottom center 3 j-symbol value.
        gamma: Bottom right 3 j-symbol value.

    Returns:
        Wigner 3 j-symbol value.
    """
    term1 = power(-1., a - b - gamma)
    numerator = sorted([x for x in [a + b - c, a - b + c, -a + b + c, a + alpha,
                                    a - alpha, b + beta, b - beta, c + gamma,
                                    c - gamma] if x > 1])
    denominator1 = [a + b + c + 1]
    v = min(a + alpha, a - alpha, b + beta, b - beta, c + gamma, c - gamma,
            a + b - c, b + c - a, c + a - b)
    s = 0.
    for t in range(v + 1):
        d = [x for x in [t, c - b + t + alpha, c - a + t - beta,
                         a + b - c - t, a - t - alpha,
                         b - t + beta] if x >= 0]
        if len(d) < 6:
            continue
        denominator = sorted([x for x in denominator1 + d + d if x > 1])
        top, bottom = [], []
        for x, y in zip(numerator, denominator):
            if x == y: continue
            if x > y:
                top.append((y + 1, x))
            else:
                bottom.append((x + 1, y))
        if len(numerator) > len(denominator):
            top += [(1, x) for x in numerator[len(denominator):]]
        elif len(numerator) < len(denominator):
            bottom += [(1, x) for x in denominator[len(numerator):]]

        f1 = 1
        for x in top:
            f1 *= _fact(x[0], x[1])
        f2 = 1
        for x in bottom:
            f2 *= _fact(x[0], x[1])
        s += power(-1., t)*sqrt(f1/f2)
    return term1*s


def _fact(start, stop):
    """Calculates the product of all integers in range start:stop (inclusive).

    Args:
        start: Range lower bound.
        stop: Range upper bound.
    """
    x = start
    for i in range(start + 1, stop + 1):
        x *= i
    return x
