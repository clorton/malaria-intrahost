#!/usr/bin/python

import duration_distribution as dist
import identity
import probability
import random


class InfectionStateChange:
    NONE = 0,
    Cleared = 1,
    Fatal = 2,
    New = 3


class MortalityTimeCourse(object):
    DAILY_MORTALITY = 0
    MORTALITY_AFTER_INFECTIOUS = 1


class Config(object):

    vital_disease_mortality = False
    incubation_distribution = dist.DurationDistribution(probability.DistributionFunction.FIXED_DURATION)
    infectious_distribution = dist.DurationDistribution(probability.DistributionFunction.FIXED_DURATION)
    base_infectivity = 1.0
    base_mortality = 1.0
    mortality_time_course = MortalityTimeCourse.DAILY_MORTALITY

    @classmethod
    def Configure(cls, values):
        for key, value in values.items():
            setattr(cls, key, value)
        return


class Infection(object):
    def __init__(self, parent=None):
        self.parent = parent
        self.duration = 0.0
        self.total_duration = 0.0
        self.incubation_timer = 0.0
        self.infectious_timer = 0.0
        self.infectiousness = 0.0
        self.StateChange = InfectionStateChange.NONE
        self.infection_strain = None
        return

    def CreateInfectionStrain(self, infstrain):
        if self.infection_strain is None:
            self.infection_strain = identity.StrainIdentity()
        if infstrain is not None:
            self.infection_strain.antigen_id = infstrain.antigen_id
            self.infection_strain.sub_strain = infstrain.sub_strain
        return

    def EvolveStrain(self, immunity, dt):
        return

    def GetInfectiousness(self):
        return self.infectiousness

    def InitInfectionImmunology(self, immunity):
        return

    def SetParameters(self, infstrain=None, incubation_period_override=-1):
        self.CreateInfectionStrain(infstrain)
        if incubation_period_override != -1:
            self.incubation_timer = incubation_period_override
        else:
            self.incubation_timer = Config.incubation_distribution.CalculateDuration()
        self.infectious_timer = Config.infectious_distribution.CalculateDuration()
        self.total_duration = self.incubation_timer + self.infectious_timer
        self.infectiousness = 0.0
        self.StateChange = InfectionStateChange.NONE
        if self.incubation_timer <= 0:
            self.infectiousness = Config.base_infectivity
        return

    def Update(self, dt, immunity):
        self.StateChange = InfectionStateChange.NONE
        self.duration += dt
        if self.duration > self.incubation_timer:
            self.infectiousness = Config.base_infectivity
        idvie = None
        vdm = Config.vital_disease_mortality
        if vdm and (Config.mortality_time_course == MortalityTimeCourse.DAILY_MORTALITY) and (self.duration > self.incubation_timer):
            idvie = self.parent.GetInterventionsContext().getDrugVaccineInterventionEffects()
            if random.random() < (Config.base_mortality * dt * immunity.getModMortality() * idvie.getInterventionReducedMortality()):
                self.StateChange = InfectionStateChange.Fatal
        if self.duration > self.total_duration:
            if vdm and (Config.mortality_time_course == MortalityTimeCourse.MORTALITY_AFTER_INFECTIOUS):
                idvie = self.parent.GetInterventionsContext().getDrugVaccineInterventionEffects()
                if random.random() < (Config.base_mortality * immunity.getModMortality() * idvie.GetInterventionReducedMortality()):
                    self.StateChange = InfectionStateChange.Fatal
                else:
                    self.StateChange = InfectionStateChange.Cleared
            else:
                self.StateChange = InfectionStateChange.Cleared
        self.EvolveStrain(immunity, dt)
        return
