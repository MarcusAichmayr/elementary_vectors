from vectors_in_intervals import *
from numpy import argmin
m_A = 6
m_B = 6
m_C = 2
length = 8
ring = ZZ

fmt = "{:10.8f}"

# for i in range(100):
#     if i % 10 == 0:
#         print(i)
#     A = random_matrix(ring, m_A, length)
#     B = random_matrix(ring, m_B, length)
#     C = random_matrix(ring, m_C, length)

#     S0 = new.AlternativesHomogeneous(A, B, C)
#     S1 = new.AlternativesHomogeneous(A, B, C)
#     S2 = old.HomogeneousSystem(A, B, C)

#     r0 = S1.certify(random=False)
#     r1 = S1.certify(random=True)
#     r2 = S2.certify()

#     if r0[0] != r1[0]:
#         raise ValueError
#     if r1[0] != r2[0]:
#         raise ValueError

#     if True:
#         print(
#             r0[0],
#             "{:>5}".format(r0[2]),
#             "{:>5}".format(r1[2]),
#             sep="; "
#         )
#     else:
#         t0 = timeit("S0.certify(random=False)", seconds=True)
#         t1 = timeit("S1.certify(random=True)", seconds=True)
#         t2 = timeit("S2.certify()", seconds=True)

#         print(
#             r1[0],
#             fmt.format(t0),
#             fmt.format(t1),
#             fmt.format(t2),
#             # fmt.format(t0 / t1),
#             # fmt.format(t1 / t2),
#             fmt.format(t0 / t2),
#             sep="; "
#         )


for i in range(100):
    A = random_matrix(ring, m_A, length)
    B = random_matrix(ring, m_B, length)
    C = random_matrix(ring, m_C, length)
    S = AlternativesHomogeneous(A, B, C)

    results = [
        S.certify(),
        S.certify(reverse=False),
        S.certify(random=True),
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
