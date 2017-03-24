#!/usr/bin/python

import identity
import infection_malaria


class Config(object):

    superinfection = False
    max_ind_inf = 1

    @classmethod
    def Update(cls, values):
        for key, value in values.items():
            setattr(cls, key, value)

        return


class Individual(object):
    def __init__(self):

        self.cumulativeInfs = 0
        self.m_is_infected = False
        self.susceptibility = None
        self.infections = []
        self.infectiousness = 0.0

        return

    def AcquireNewInfection(self, cp, incubation_period_override):
        newStrainId = identity.StrainIdentity()
        if cp is not None:
            newStrainId = cp.ResolveInfectingStrain()

        numInfs = len(self.infections)
        if (numInfs == 0) or (Config.superinfection and (numInfs < Config.max_ind_inf)):
            self.cumulativeInfs += 1
            self.m_is_infected = True

            newInf = self.createInfection()
            newInf.SetParameters(newStrainId, incubation_period_override)
            newInf.InitInfectionImmunology(self.susceptibility)

            self.infections.append(newInf)
            self.infectiousness += newInf.GetInfectiousness()

            # self.ReportInfectionState()

        return

    def createInfection(self):
        return infection_malaria.InfectionMalaria()

    def GetInterventionsContext(self):
        return None

    def Update(self,currenttime, dt):
        return