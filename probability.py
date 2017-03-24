#!/usr/bin/python

import math
import random
import randgen


class DistributionFunction(object):
    NOT_INITIALIZED = -1
    FIXED_DURATION = 0
    UNIFORM_DURATION = 1
    GAUSSIAN_DURATION = 2
    EXPONENTIAL_DURATION = 3
    POISSON_DURATION = 4
    LOG_NORMAL_DURATION = 5
    BIMODAL_DURATION = 6
    PIECEWISE_CONSTANT = 7
    PIECEWISE_LINEAR = 8
    WEIBULL_DURATION = 9


def fromDistribution(function, param1, param2, default_value):
    value = 0.0
    if function == DistributionFunction.FIXED_DURATION:
        value = param1
    elif function == DistributionFunction.UNIFORM_DURATION:
        value = param1 + random.random() * (param2 - param1)
    elif function == DistributionFunction.GAUSSIAN_DURATION:
        value = randgen.eGauss()* param2 + param1
        if value < 0:
            value = 0
    elif function == DistributionFunction.EXPONENTIAL_DURATION:
        value = randgen.expdist(param1)
    elif function == DistributionFunction.POISSON_DURATION:
        value = randgen.Poisson(param1)
    elif function == DistributionFunction.LOG_NORMAL_DURATION:
        value = math.exp(math.log(param1) + randgen.eGauss() * param2)
    elif function == DistributionFunction.BIMODAL_DURATION:
        value = param2 if random.random() < param1 else 1.0
    elif function == DistributionFunction.WEIBULL_DURATION:
        value = randgen.Weibulil(param1, param2)
    else:
        assert(False)
    return value