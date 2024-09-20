from vectors_in_intervals import *
from numpy import argmin
m_A = 6
m_B = 6
length = 5
ring = ZZ

fmt = "{:10.8f}"

print(
    "{:<2}".format("i"),
    "{:<5}".format("res"),
    *[
        "{:<5}".format(n) for n in [
            "s0",
            "rev0",
            "rand0",
            "s1",
            "rev1",
            "rand1",
            "s2",
            "rev2",
            "rand2",
        ]
    ],
    sep="; "
)

for i in range(100):
    A = random_matrix(ring, m_A, length)
    B = random_matrix(ring, m_B, length)
    b = random_vector(ring, m_A)
    c = random_vector(ring, m_B)
    S0 = AlternativesInhomogeneous(A, B, b, c)
    S1 = AlternativesInhomogeneous(A, B, b, c, one_homogenized=True)
    S2 = AlternativesInhomogeneous(A, B, b, c, two_double_system=True)

    results = [
        S0.certify(),
        S0.certify(reverse=False),
        S0.certify(random=True),
        S1.certify(),
        S1.certify(reverse=False),
        S1.certify(random=True),
        S2.certify(),
        S2.certify(reverse=False),
        S2.certify(random=True),
    ]

    if all(r[0] for r in results) != any(r[0] for r in results):
        raise ValueError("Wrong result!")

    if True:
        print(
            "{:>2}".format(i),
            "{:>5}".format(str(results[0][0])),
            *["{:>5}".format(r[2]) for r in results],
            argmin([r[2] for r in results]),
            sep="; "
        )
    else:
        kwargs = {"seconds": True, "number": 1, "repeat": 1}
        times = [
            # timeit("S0.certify()", **kwargs),
            timeit("S1.certify()", **kwargs),
            # timeit("S2.certify()", **kwargs),
            timeit("S3.certify()", **kwargs),
        ]

        print(
            "{:>5}".format(str(results[0][0])),
            *[fmt.format(ti) for ti in times],
            *["{:>5}".format(r[2]) for r in results],
            argmin(times),
            sep="; "
        )
