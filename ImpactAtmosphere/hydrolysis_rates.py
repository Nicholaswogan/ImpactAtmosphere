

def HCN_hydrolysis_rate(T, pH):
    '''Calculates HCN hydrolysis rate following Miyakawa et al. 2001.

    Parameters
    ----------
    T: float
        Temperature of the water (K)
    pH; float
        The pH of the water (unit-less)

    Returns
    -------
    ktot: float
        The HCN hydrolysis rate (1/s)
    '''
    H=10.**(-1.*pH)
    OH=1.0e-14/H

    # acid catalyzed: (in M^-1 s^-1)
    logk1H=-4950./T + 8.43
    k1H=10.**(logk1H)
    # base catalyzed: (in M^-1 s^-1)
    logk1OH=-4240./T + 11.1
    k1OH=10.**(logk1OH)

    pKw=-6.0846+4471.33/T+.017053*T  #from Stribling and Miller 1987, from Miyakawa's original sources (Robinson and Stokes 1959 and Schlesinger and Miller 1973, resp.)
    pKa_HCN=-8.85+3802./T+.01786*T
    Kw=10.**(-1.*pKw)
    Ka_HCN=10.**(-1.*pKa_HCN)

    ktot=k1H*H+(k1OH*Kw/(H+Ka_HCN))  #in s^-1

    return ktot # in s^-1

def HCONH2_hydrolysis_rate(T,pH):
    '''Calculates HCONH2 hydrolysis rate following Miyakawa et al. 2001.

    Parameters
    ----------
    T: float
        Temperature of the water (K)
    pH; float
        The pH of the water (unit-less)

    Returns
    -------
    ktot: float
        The HCONH2 hydrolysis rate (1/s)
    '''
    H=10.**(-1.*pH)
    OH=1.0e-14/H

    logkformH=-4060./T+9.85
    logkformOH=-3440./T+8.92

    kformH=10.**(logkformH)
    kformOH=10.**(logkformOH)

    kform=kformH*H+kformOH*OH

    return kform
