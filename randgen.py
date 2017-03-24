#!/usr/bin/python

import math
import random


def eGauss():
    rad = -2 * math.log(random.random())

    r1 = random.random() - 0.5
    r2 = random.random() - 0.5
    s = r1 * r1 + r2 * r2

    while s > 0.25:

        r1 = random.random() - 0.5
        r2 = random.random() - 0.5
        s = r1*r1 + r2*r2

    norm = math.sqrt(rad/s)

    return r2 * norm


def expdist(rate):
    return -math.log(random.random()) / rate if rate != 0 else 0.0


def Poisson(ratetime):
    if ratetime <= 0:
        return 0.0

    events = 0

    if ratetime < 10:
        Time = 0
        while Time < 1:
            Time += math.log(random.random()) / ratetime
            if Time < 1:
                events += 1
    else:
        tempval = eGauss() * math.sqrt(ratetime) + ratetime + 0.5
        events = int(tempval) if tempval >= 0 else 0

    return events


def Weibull(lambduh, kappa):
    return lambduh * math.pow( -math.log(random.random()), 1.0/kappa)