from mpmath import mp
from fresnel import fresnelConstalesLog, log_ball_vol


def volumeHypersphericalCap(N, h):
    """
    Using (3) Concise Formulas for the Area and Volume of a Hyperspherical Cap
    S. Li, Asian Journal of Mathematics & Statistics, 4: 66-70.

    :param N: we consider a hypersphere in RR^N
    :param h: the height of the spherical cap, i.e. perpendicular distance from
                the plane that intersects hypersphere to surface of hypersphere

    :returns: the volume of a hyperspherical cap of height h

    ..note:: the article works in terms of phi, which is related to h via
                h = (1 - cos phi) r
    ..note:: we are assuming r = 1 throughout
    """
    assert N >= 2
    assert 0 <= h <= 1

    mpN = mp.mpf(N)
    mph = mp.mpf(h)
    mp1 = mp.mpf(1)
    mp2 = mp.mpf(2)
    mp3 = mp.mpf(3)

    phi = mp.acos(mp.mpf(1) - mph)

    # I_(sin^2 phi)((N + 1)/2, 1/2)
    regularisedParam = mp.power(mp.sin(phi), mp2)
    betaTerm = mp.betainc((mpN + mp1) / mp2, mp1 / mp2, x2=regularisedParam,
                          regularized=True)
    capVolume = mp.power(mp.e, log_ball_vol(mpN)) * betaTerm / mp2

    # sanity checks
    if N == 2:
        correctValue = phi - mp.sin(phi)*mp.cos(phi)
        assert abs(correctValue - capVolume) <= 10**-12
    if N == 3:
        correctValue = (mp2/mp3 - mp.cos(phi) + mp.power(mp.cos(phi), mp3)/mp3) * mp.pi # noqa
        assert abs(correctValue - capVolume) <= 10**-12

    return capVolume


def testSmalls(N, s=1.5):
    """
    For N >=2 and s in (1, 2) only the dimension N - 1 spherical caps are not
    contained within the cube. There are 2N such and volumeHypersphericalCap
    gives us the volume of one such.

    In this instance, we can check fresnelConstalesLog against the ground truth
    for growing dimensions in the case where it appears to become inaccurate
    most quickly, i.e. small s. This (heuristically) acts as a worst case.

    :param N: we work in RR^N
    :param s: a scaling parameter, here taken in (1, 2)
    """
    assert 1 < s < 2

    mpN = mp.mpf(N)
    mps = mp.mpf(s)
    mp1 = mp.mpf(1)
    mp2 = mp.mpf(2)

    h = mp1 - mp1 / mp.sqrt(mps)
    capsVolume = mp2 * mpN * volumeHypersphericalCap(N, h)
    intersectionVolume = mp.power(mp.e, log_ball_vol(N)) - capsVolume

    fresnelVolume = mp.power(mp.e, fresnelConstalesLog(N, s))

    assert abs(intersectionVolume - fresnelVolume) <= 10**-5

    return intersectionVolume, fresnelVolume


def findAccurateRange():
    accurateN = None
    worsts = None
    ss = [1 + 0.05 * i for i in range(1, 20)]
    for s in ss:
        for N in range(200, 351, 5):
            try:
                testSmalls(N, s=s)
            # the exact answer differs too much from the fresnel estimate
            except AssertionError:
                previousN = N - 5
                if accurateN is None or previousN < accurateN:
                    accurateN = previousN
                    worsts = s
                break
            print("good (s, N):", s, N)
    return accurateN, worsts
