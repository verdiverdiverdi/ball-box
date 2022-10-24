import os
import pickle

from mpmath import mp

default_prec = 12 * 53
mp.prec = default_prec


def default_terms(N):
    return int(max(10*N, default_prec, 10000))


def log_ball_vol(d, R=mp.mpf(1), prec=None):
    if prec is None:
        mp.prec = default_prec
    mpd = mp.mpf(d)
    return mpd * mp.log(R) + (mpd / mp.mpf(2)) * mp.log(mp.pi) - mp.loggamma(mpd/mp.mpf(2) + mp.mpf(1)) # noqa


def fresnelConstalesLog(N, s, terms=None, prec=None):
    """
    Following (2.6) of https://arxiv.org/abs/1804.07861
    Comment on "Sum of squares of uniform random variables" by I. Weissman,
    Peter J. Forrester, we calculate an estimate of F_N(s) via truncated
    infinite series.

    This is an estimate for vol([0, 1]^N ∩ B_N(sqrt(s))), so that
    2^N s^(-N/2) F_N(s) = vol([-1/sqrt(s), 1/sqrt(s)]^N ∩ B_N(1)).

    Therefore, if R is the radius of our ball and q the side length of our cube
    we set s such that q/(2R) = 1/sqrt(s) in our use cases.

    :param N:       we work in RR^N
    :param s:       a scaling constant
    :param terms:   the number of terms in the infinite sum to consider,
                        if ``None`` use default_terms(N)
    :param prec:    the amount of precision with which to compute with, see
                        mpmath documentation

    :returns:       log(vol([-1/sqrt(s), 1/sqrt(s)]^N ∩ B_N(1)))
    """
    if prec is None:
        mp.prec = default_prec
    else:
        try:
            mp.prec = prec
        except TypeError:
            raise Exception("Incorrect type for precision.")

    mpN = mp.mpf(N)
    mps = mp.mpf(s)

    if mps <= mp.mpf('1'):
        # the entire ball
        return log_ball_vol(mpN)

    if mps >= mpN:
        # the entire cube
        return N * (mp.log(mp.mpf('2')) - mp.mpf(mp.mpf('1') / mp.mpf('2')) * mp.log(mps))  # noqa

    constants = mp.mpf('1') / mp.mpf('6') + mps / mpN

    if terms is None:
        terms = default_terms(N)

    fullSum = mp.mpf('0')

    precomputeFresnelTerms(N, prec=mp.prec, terms=terms)
    dimNPrecomputed = retrievePrecomputed(N, prec=mp.prec, terms=terms)

    for k in range(1, terms + 1):
        mpk = mp.mpf(str(int(k)))

        fresnelFinal = dimNPrecomputed[k]

        if mpk * mps / mpN % mp.mpf('1') == mp.mpf('0'):
            eTerm = mp.mpf('1') / mpk
        else:
            exponent = mp.mpf('2') * mp.pi * mp.j * mpk * mps / mpN
            eTermPre = mp.power(mp.e, exponent)
            eTerm = eTermPre / mpk

        summand = fresnelFinal * eTerm
        fullSum += summand

    imaginarySum = fullSum.imag
    imaginarySumFinal = imaginarySum / mp.pi
    # we have calculated (2.6)
    FNs = constants + imaginarySumFinal

    assert FNs > mp.mpf('0'), "either more terms or precision required"
    # we return log(vol([-1/sqrt(s), 1/sqrt(s)]^N ∩ B_N(1))) using (2.3)
    correctionFactor = mpN * mp.log(mp.mpf('2') / mp.sqrt(mps))
    return correctionFactor + mp.log(FNs)


def retrievePrecomputed(N, terms=None, prec=None):
    if terms is None:
        terms = default_terms(N)
    if prec is None:
        prec = default_prec

    filename = "precomps/" + "dim{N}-prec{p}".format(N=str(N), p=str(prec))

    try:
        with open(filename, 'rb') as fh:
            precomputed = pickle.load(fh)
    except FileNotFoundError:
        raise ("Ensure precomputation for dimension {N} and precision {p}"
               .format(N=str(N), p=str(prec)))

    # check all required ks are present in dimNPrecomputed
    if not all(x in precomputed.keys() for x in range(1, terms + 1)):
        raise ("Some k are missing for dimension {N} and precision {p}"
               .format(N=str(N), p=str(prec)))

    return precomputed


def precomputeFresnelTerms(N, terms=None, prec=None, verbose=False):
    """
    Precomputes (per k in (2.6) of https://arxiv.org/abs/1804.07861) the term
        - (1 / k) * (C(x) - i S(x) / (x))**N
    where x = sqrt(2 * pi * k / N), and saves them in high precision as a
    dictionary

    :param N:       we work in RR^N
    :param terms:   the number of terms in the infinite sum to consider,
                        if ``None`` use default_terms(N)
    :param prec:    the amount of precision with which to compute with, see
                        mpmath documentation

    ..note::        this is assuming the fixed version of the arXiv paper is
                    uploaded

    ..note::        there are two definitions of the Fresnel integrals used,
                    the first is int_0^x (sin(t^2) dt), as in (2.6), and
                    similarly for cos, then there is the ``normalised``
                    int_0^x (sin(pi t^2 / 2) dt), and similarly for cos.

                    If we want C(x) and S(x) then:
                        - C(x) = sqrt(pi / 2) * mp.fresnelc(sqrt(2 / pi) * x)
                        - S(x) = sqrt(pi / 2) * mp.fresnels(sqrt(2 / pi) * x)
                    hence the fiddling that goes on below
    """
    if terms is None:
        terms = default_terms(N)
    if prec is None:
        mp.prec = default_prec
    else:
        mp.prec = prec

    if not os.path.exists("precomps/"):
        os.mkdir("precomps/")

    filename = "precomps/" + "dim{N}-prec{p}".format(N=str(N), p=str(mp.prec))

    try:
        with open(filename, 'rb') as fh:
            precomputed = pickle.load(fh)
    except FileNotFoundError:
        if verbose:
            print("Creating:", filename)
        precomputed = {}

    # test whether all necessary terms already precomputed
    if all(x in precomputed.keys() for x in range(1, terms+1)):
        return

    mpN = mp.mpf(N)

    for k in range(1, terms + 1):
        if k in precomputed.keys():
            continue

        mpk = mp.mpf(str(int(k)))

        fresnelInput = mp.sqrt(mp.mpf('2') * mp.pi * mpk / mpN)
        fresnelCorrection = mp.sqrt(mp.mpf('2') / mp.pi)
        fresnelcPre = mp.fresnelc(fresnelCorrection * fresnelInput)
        fresnelsPre = mp.fresnels(fresnelCorrection * fresnelInput) * mp.j
        fresnelPreEstimate = mp.sqrt(mp.pi / mp.mpf('2')) * (fresnelcPre - fresnelsPre)  # noqa
        fresnelEstimate = fresnelPreEstimate / fresnelInput
        fresnelFinal = mp.power(fresnelEstimate, mpN)

        precomputed[k] = fresnelFinal

    with open(filename, 'wb') as fh:
        pickle.dump(precomputed, fh, protocol=pickle.HIGHEST_PROTOCOL)


def asymptoticEstimate(N, s):
    """
    From Thm 3 of https://eprint.iacr.org/2017/155.pdf, for balanced boxes.
    Gives an asymptotic estimate for logvol(B_N(1) ∩ [-1/sqrt(s),1/sqrt(s)]^N)

    ..note:: method only works for us when s > N/3
    """
    assert s > N/3
    mpN = mp.mpf(N)
    mps = mp.mpf(s)
    # first term of (3)
    logBox = mpN * mp.log(mp.mpf(2) / mp.sqrt(mps))
    # appealing to Thm 3 itself
    Exp = mpN / (mp.mpf(3) * mps)
    Var = mp.mpf(4) * mpN / (mp.mpf(45) * mp.power(s, mp.mpf(2)))
    # solve for y
    mpy = (mp.mpf(1) - Exp) / Var
    prob = (mp.mpf(1) + mp.erf(mpy / mp.sqrt(mp.mpf(2)))) / mp.mpf(2)
    return logBox + mp.log(prob)
