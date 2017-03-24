#!/usr/bin/python

import probability as Probability


class DurationDistribution(object):
    def __init__(self, default_type=Probability.DistributionFunction.NOT_INITIALIZED):

        self.m_Type = default_type
        self.m_Param1 = 0.0
        self.m_Param2 = 0.0

        return

    def CalculateDuration(self):
        p1 = self.m_Param1
        p2 = self.m_Param2

        if self.m_Type == Probability.DistributionFunction.EXPONENTIAL_DURATION:
            p1 = 1.0/p1
        elif self.m_Type == Probability.DistributionFunction.WEIBULL_DURATION:
            p2 = 1.0/p2

        duration = Probability.fromDistribution(self.m_Type, p1, p2, 0.0)

        return duration
