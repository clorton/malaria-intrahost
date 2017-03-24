#!/usr/bin/python

import susceptibility_malaria
import sigmoid


NON_TRIVIAL_ANTIBODY_THRESHOLD = 0.0000001
TWENTY_DAY_DECAY_CONSTANT = 0.05
B_CELL_PROLIFERATION_THRESHOLD = 0.4
B_CELL_PROLIFERATION_CONSTANT = 0.33
ANTIBODY_RELEASE_THRESHOLD = 0.3
ANTIBODY_RELEASE_FACTOR = 4


class AntibodyType:
    CSP = 0
    MSP1 = 1
    PfEMP1_minor = 2
    PfEMP1_major = 3
    N_MALARIA_ANTIBODY_TYPES = 4


class Antibody(object):
    def __init__(self, antibody_type, variant, capacity=0, concentration=0):

        self.m_antigen_count = 0
        self.m_antigen_present = False

        self.m_antibody_type = antibody_type
        self.m_antibody_variant = variant
        self.m_antibody_capacity = capacity
        self.m_antibody_concentration = concentration

        return

    def Decay(self, dt):

        if self.m_antibody_concentration > NON_TRIVIAL_ANTIBODY_THRESHOLD:
            self.m_antibody_concentration -= self.m_antibody_concentration * TWENTY_DAY_DECAY_CONSTANT * dt

        if self.m_antibody_capacity > susceptibility_malaria.Config.memory_level:
            self.m_antibody_capacity -= (self.m_antibody_capacity - susceptibility_malaria.Config.memory_level) * susceptibility_malaria.Config.hyperimmune_decay_rate * dt

        return

    def StimulateCytokines(self, dt, inv_uL_blood):

        return (1 - self.m_antibody_concentration) * self.m_antigen_count * inv_uL_blood

    def UpdateAntibodyCapacity(self, dt, inv_uL_blood):

        growth_rate = susceptibility_malaria.Config.MSP1_antibody_growthrate
        threshold = susceptibility_malaria.Config.antibody_stimulation_c50

        self.m_antibody_capacity += growth_rate * (1.0 - self.m_antibody_capacity) * sigmoid.basic_sigmoid(threshold, self.m_antigen_count * inv_uL_blood)

        if self.m_antibody_capacity > B_CELL_PROLIFERATION_THRESHOLD:
            self.m_antibody_capacity += (1.0 - self.m_antibody_capacity) * B_CELL_PROLIFERATION_CONSTANT * dt

        if self.m_antibody_capacity > 1:
            self.m_antibody_capacity = 1.0

        return

    def UpdateAntibodyCapacityByRate(self, dt, growth_rate):

        self.m_antibody_capacity += growth_rate * dt * (1.0 - self.m_antibody_capacity)

        if self.m_antibody_capacity > 1:
            self.m_antibody_capacity = 1.0

        return

    def UpdateAntibodyConcentration(self, dt):

        if self.m_antibody_capacity > ANTIBODY_RELEASE_THRESHOLD:
            self.m_antibody_concentration += (self.m_antibody_capacity - self.m_antibody_concentration) * ANTIBODY_RELEASE_FACTOR * dt

        if self.m_antibody_concentration > self.m_antibody_capacity:
            self.m_antibody_concentration = self.m_antibody_capacity

        return

    def ResetCounters(self):

        self.m_antigen_present = False
        self.m_antigen_count = 0

        return

    def IncreateAntigenCount(self, antigenCount):

        if antigenCount > 0:
            self.m_antigen_count += antigenCount
            self.m_antigen_present = True

        return


class AntibodyCSP(Antibody):
    def __init__(self, variant, capacity):
        super(AntibodyCSP, self).__init__(AntibodyType.CSP, variant, capacity)
        return

    def Decay(self, dt):

        if self.m_antibody_concentration > self.m_antibody_capacity:
            self.m_antibody_concentration -= self.m_antibody_capacity * dt / susceptibility_malaria.Config.antibody_csp_decay_days
        else:
            super(AntibodyCSP, self).Decay(dt)

        return

    def UpdateAntibodyConcentration(self, dt):

        if self.m_antibody_concentration > self.m_antibody_capacity:
            self.m_antibody_concentration -= self.m_antibody_concentration * dt / susceptibility_malaria.Config.antibody_csp_decay_days
        else:
            super(AntibodyCSP, self).UpdateAntibodyConcentration(dt)

        return


class AntibodyMSP(Antibody):
    def __init__(self, variant, capacity):
        super(AntibodyMSP, self).__init__(AntibodyType.MSP1, variant, capacity)
        return


class AntibodyPfEMP1Minor(Antibody):
    def __init__(self, variant, capacity):
        super(AntibodyPfEMP1Minor, self).__init__(AntibodyType.PfEMP1_minor, variant, capacity)
        return

    def UpdateAntibodyCapacity(self, dt, inv_uL_blood):

        min_stimulation = susceptibility_malaria.Config.antibody_stimulation_c50 * susceptibility_malaria.Config.minimum_adapted_response
        growth_rate = susceptibility_malaria.Config.antibody_capacity_growthrate * susceptibility_malaria.Config.non_specific_growth
        threshold = susceptibility_malaria.Config.antibody_stimulation_c50

        if self.m_antibody_capacity <= B_CELL_PROLIFERATION_THRESHOLD:
            self.m_antibody_capacity += growth_rate * dt * (1.0 - self.m_antibody_capacity) * sigmoid.basic_sigmoid(threshold, self.m_antigen_count * inv_uL_blood * min_stimulation)
        else:
            self.m_antibody_capacity += (1.0 - self.m_antibody_capacity) * B_CELL_PROLIFERATION_CONSTANT * dt

        if self.m_antibody_capacity > 1:
            self.m_antibody_capacity = 1.0

        return


class AntibodyPfEMP1Major(Antibody):
    def __init__(self, variant, capacity):
        super(AntibodyPfEMP1Major, self).__init__(AntibodyType.PfEMP1_major, variant, capacity)
        return

    def UpdateAntibodyCapacity(self, dt, inv_uL_blood):

        min_stimulation = susceptibility_malaria.Config.antibody_stimulation_c50 * susceptibility_malaria.Config.minimum_adapted_response
        growth_rate = susceptibility_malaria.Config.antibody_capacity_growthrate
        threshold = susceptibility_malaria.Config.antibody_stimulation_c50

        if self.m_antibody_capacity <= B_CELL_PROLIFERATION_THRESHOLD:
            self.m_antibody_capacity += growth_rate * dt * (1.0 - self.m_antibody_capacity) * sigmoid.basic_sigmoid(threshold, self.m_antigen_count * inv_uL_blood * min_stimulation)

            if self.m_antibody_capacity > 1:
                self.m_antibody_capacity = 1.0

        else:
            self.m_antibody_capacity += (1.0 - self.m_antibody_capacity) * B_CELL_PROLIFERATION_CONSTANT * dt

        return


class pfemp1_antibody_t(object):
    def __init__(self):
        self.minor = None
        self.major = None

        return
