#!/usr/bin/python

import math

DEFAULT_ANTIBODY_CSP_DECAY_DAYS = 90


class Config(object):

    memory_level = 0.2
    MSP1_antibody_growthrate = 0.02
    antibody_stimulation_c50 = 10.0
    antibody_capacity_growthrate = 0.1
    minimum_adapted_response = 0.02
    non_specific_growth = 0.5
    antibody_csp_decay_days = DEFAULT_ANTIBODY_CSP_DECAY_DAYS

    hyperimmune_decay_rate = -math.log((0.4 - memory_level) / (1.0 - memory_level)) / 120.0

    @classmethod
    def Update(cls, values):
        for key, value in values.items():
            setattr(cls, key, value)

        cls.hyperimmune_decay_rate = -math.log((0.4 - cls.memory_level) / (1.0 - cls.memory_level)) / 120.0

        return
