from vectors_in_intervals import *
from vectors_in_intervals.certifying_inequalities import *
from numpy import argmin
length_m = 4
length_n = 5
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
    M = random_matrix(ring, length_m, length_n)
    I = random_intervals(length_m, ring=ring)
    S0 = AlternativesGeneral(M, I)
    S1 = AlternativesInhomogeneous(*general_to_inhomogeneous(M, I))
    S2 = AlternativesHomogeneous(*inhomogeneous_to_homogeneous(*general_to_inhomogeneous(M, I)))

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
