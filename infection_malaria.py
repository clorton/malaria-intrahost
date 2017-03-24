#!/usr/bin/python

import math
import random

import antibody
import identity
import sigmoid
import randgen


class ParasiteSwitchType:
    CONSTANT_SWITCH_RATE_2VARS = 0
    RATE_PER_PARASITE_7VARS = 1
    RATE_PER_PARASITE_5VARS_DECAYING = 2


class MalariaStrains:
    FALCIPARUM_NONRANDOM_STRAIN = 11
    FALCIPARUM_RANDOM50_STRAIN = 12
    FALCIPARUM_RANDOM_STRAIN = 13
    FALCIPARUM_STRAIN_GENERATOR = 20


class Config(object):

    parasite_switch_type = ParasiteSwitchType.CONSTANT_SWITCH_RATE_2VARS
    malaria_strains = MalariaStrains.FALCIPARUM_NONRANDOM_STRAIN
    antibody_IRBC_killrate = 2.0
    MSP1_merozoite_kill = 0.5
    gametocyte_stage_survival = 1.0
    base_gametocyte_sexratio = 0.2
    base_gametocyte_production = 0.02
    antigen_switch_rate = 0.33
    merozoites_per_hepatocyte = 15000
    merozoites_per_schizont = 16
    non_specific_antigenicity = 0.2
    RBC_destruction_multiplier = 9.5
    n_asexual_cycles_wo_gametocytes = 1

    @classmethod
    def Update(cls, values):
        for key, value in values.items():
            setattr(cls, key, value)

        return


class AsexualCycleStatus:
    NoAsexualCycle = 0
    AsexualCycle = 1
    HepatocyteRelease = 2


class GametocyteStages:
    Stage0 = 0
    Stage1 = 1
    Stage2 = 2
    Stage3 = 3
    Stage4 = 4
    Mature = 5
    Count = 6


CLONAL_PfEMP1_VARIANTS = 50
INV_MICROLITERS_BLOOD_ADULT = 1.0/5e6
MEROZOITE_LIMITING_RBC_THRESHOLD = 0.2
MIN_FEVER_DEGREES_KILLING = 1.5
MINIMUM_IRBC_COUNT = 10


class InfectionMalaria(object):
    def __init__(self, initial_hepatocytes=1):

        self.m_IRBCtimer = 0.0
        self.m_hepatocytes = 0
        self.m_asexual_phase = AsexualCycleStatus.NoAsexualCycle
        self.m_asexual_cycle_count = 0
        self.m_MSPtype = 0
        self.m_nonspectype = 0
        self.m_minor_epitope_type = [0 for x in range(CLONAL_PfEMP1_VARIANTS)]
        self.m_IRBCtype = [0 for x in range(CLONAL_PfEMP1_VARIANTS)]
        self.m_MSP_antibody = None
        self.m_PfEMP1_antibodies = [antibody.pfemp1_antibody_t() for x in range(CLONAL_PfEMP1_VARIANTS)]
        self.m_IRBC_count = [0 for x in range(CLONAL_PfEMP1_VARIANTS)]
        self.m_malegametocytes = [0 for x in range(GametocyteStages.Count)]
        self.m_femalegametocytes = [0 for x in range(GametocyteStages.Count)]
        self.m_gametorate = 0.0
        self.m_gametosexratio = 0.0
        self.m_measured_duration = 0
        self.m_start_measuring = False
        self.m_temp_duration = 0
        self.m_max_parasites = 0
        self.m_inv_microliters_blood = INV_MICROLITERS_BLOOD_ADULT
        self.drugResistanceFlag = 0
        self.m_pMDE = None

        self.duration = 0
        self.infectiousness = 0
        self.StateChange = InfectionStateChange.NONE

        self.m_hepatocytes = initial_hepatocytes
        
        self.parent = None
        self.infection_strain = None

        return

    def Update(self, dt, immunity):

        self.m_inv_microliters_blood = immunity.get_inv_microliters_blood()

        self.StateChange = InfectionStateChange.NONE
        self.duration += dt

        if self.m_hepatocytes > 0:
            self.malariaProcessHepatocytes(dt, immunity)

        if self.m_asexual_phase > AsexualCycleStatus.NoAsexualCycle:
            if self.m_asexual_phase == AsexualCycleStatus.HepatocyteRelease:
                self.m_asexual_phase = AsexualCycleStatus.AsexualCycle
            else:
                self.m_IRBCtimer -= dt

            if self.m_IRBCtimer <= 0:
                self.processEndOfAsexualCycle(immunity)

            if immunity.get_RBC_count() < 1:
                # individual dies
                pass

            self.malariaImmuneStimulation(dt, immunity)
            self.malariaImmunityIRBCKill(dt, immunity)
            self.malariaImmunityGametocyteKill(dt, immunity)

            self.m_MSP_antibody.IncreateAntigenCount(1)
            immunity.SetAntigenPresent()

        self.malariaCheckinfectionStatus(dt, immunity)

        return

    def malariaProcessHepatocytes(self, dt, immunity):

        if (dt > 0) and (self.m_hepatocytes > 0):
            context = self.parent.GetInterventionsContext()

            if self.m_pMDE is None:
                # Get m_pMDE from context
                pass

            drug_killrate = 0
            if self.m_pMDE is not None:
                drug_killrate = self.m_pMDE.get_drug_hepatocyte()

            if drug_killrate > 0:
                tempval1 = 0
                pkill = 1 - math.exp(-dt * drug_killrate)
                for x in range(self.m_hepatocytes):
                    if random.random() < pkill:
                        tempval1 += 1

                self.m_hepatocytes -= tempval1

            incubation_period = Config.incubation_distribution.GetParam1()
            if (self.m_asexual_phase == AsexualCycleStatus.NoAsexualCycle) and (self.duration >= incubation_period):
                self.m_IRBC_count.assign(CLONAL_PfEMP1_VARIANTS, 0)

                INITIAL_PFEMP1_VARIANTS = 5
                for i in range(INITIAL_PFEMP1_VARIANTS):
                    self.m_IRBC_count[i] = self.m_hepatocytes * Config.merozoites_per_hepatocyte / INITIAL_PFEMP1_VARIANTS
                    immunity.UpdateActiveAntibody( self.m_PfEMP1_antibodies[i], self.m_minor_epitope_type[i], self.m_IRBCtype[i] )

                self.m_hepatocytes = 0
                self.m_IRBCtimer = 2.0
                self.m_asexual_phase = AsexualCycleStatus.HepatocyteRelease

        return

    def processEndOfAsexualCycle(self, immunity):

        RBCavailability = immunity.get_RBC_availability()
        merozoitesurvival = max(0.0, (1.0 - Config.MSP1_merozoite_kill * self.m_MSP_antibody.GetAntibodyConcentration()) * math.exp(-RBCavailability / MEROZOITE_LIMITING_RBC_THRESHOLD))
        totalIRBC = sum(self.m_IRBC_count)
        self.m_MSP_antibody.IncreaseAntigenCount(totalIRBC)
        self.malariaCycleGametocytes(merozoitesurvival)
        self.malariaIRBCAntigenSwitch(merozoitesurvival)
        totalIRBC = 0
        for j in range(CLONAL_PfEMP1_VARIANTS):
            if self.m_IRBC_count[j] > 0:
                totalIRBC += self.m_IRBC_count[j]
                immunity.UpdateActiveAntibody(self.m_PfEMP1_antibodies[j], self.m_minor_epitope_type[j], self.m_IRBCtype[j])

        destruction_factor = max(1.0, Config.RBC_destruction_multiplier * math.exp(-RBCavailability / MEROZOITE_LIMITING_RBC_THRESHOLD))
        immunity.remove_RBCs(totalIRBC, self.m_malegametocytes[0] + self.m_femalegametocytes[0], destruction_factor)
        self.m_IRBCtimer = 2.0
        self.m_asexual_cycle_count += 1

        return

    def malariaImmuneStimulation(self, dt, immunity):

        assert((dt > 0) and (immunity is not None))

        for i in range(CLONAL_PfEMP1_VARIANTS):
            if self.m_IRBC_count[i] < 0:
                self.m_IRBC_count[i] = 0

            if self.m_IRBC_count[i] > 0:
                self.m_PfEMP1_antibodies[i].major.IncreaseAntigenCount(self.m_IRBC_count[i])
                self.m_PfEMP1_antibodies[i].minor.IncreaseAntigenCount(self.m_IRBC_count[i])
                immunity.SetAntigenPresent()

        return

    def malariaImmunityIRBCKill(self, dt, immunity):

        if (dt <= 0) or immunity is None:
            return

        fever_cytokine_killrate = immunity.get_fever_killing_rate() * sigmoid.basic_sigmoid(1.0, immunity.get_fever() - MIN_FEVER_DEGREES_KILLING) if (immunity.get_fever() > MIN_FEVER_DEGREES_KILLING) else 0.0

        patient = self.GetParent()
        context = patient.GetInterventionsContext()

        if self.m_pMDE is None:
            self.m_pMDE = context.get_IMalariaDrugEffects()

        drug_killrate = self.m_pMDE.get_drug_IRBC_killrate() if self.m_pMDE or (self.getDrugResistanceFlag() > 0) else 0

        for i in range(CLONAL_PfEMP1_VARIANTS):
            if self.m_IRBC_count[i] > 0:
                pkill = math.exp(-dt*((self.m_PfEMP1_antibodies[i].major.GetAntibodyConcentration()+Config.non_specific_antigenicity*self.m_PfEMP1_antibodies[i].minor.GetAntibodyConcentration()+immunity.get_maternal_antibodies())*Config.antibody_IRBC_killrate+fever_cytokine_killrate+drug_killrate))

                tempval1 = self.m_IRBC_count[i]*pkill
                if tempval1 > 0:
                    tempval1 = randgen.eGauss()*math.sqrt(tempval1*(1.0-pkill))+tempval1

                if tempval1 < 0.5:
                    tempval1 = 0

                self.m_IRBC_count[i] -= int(tempval1+0.5)

                if self.m_IRBC_count[i] < 1:
                    self.m_IRBC_count[i] = 0

        return

    def malariaImmunityGametocyteKill(self, dt, immunity):

        patient = self.GetParent()
        context = patient.GetInterventionsContext()

        if (dt <= 0) or not immunity:
            return

        for i in range(GametocyteStages.Mature):
            fever_cytokine_killrate = 0
            if not self.m_pMDE:
                self.m_pMDE = context.get_IMalariaDrugEffects()
            drug_killrate = 0
            if self.m_pMDE:
                if i < GametocyteStages.Stage3:
                    drug_killrate = self.m_pMDE.get_drug_gametocyte02()
                else:
                    drug_killrate = self.m_pMDE.get_drug_gametocyte34()
            gametocyte_kill_fraction = math.exp(-dt*(fever_cytokine_killrate+drug_killrate))
            self.m_malegametocytes[i] -= int(0.5+self.m_malegametocytes[i]*gametocyte_kill_fraction)
            if self.m_malegametocytes[i] < 1:
                self.m_malegametocytes[i] = 0
            self.m_femalegametocytes[i] -= int(0.5 + self.m_femalegametocytes[i] * gametocyte_kill_fraction)
            if self.m_femalegametocytes[i] < 1:
                self.m_femalegametocytes[i] = 0

        return

    def malariaCheckinfectionStatus(self, dt, immunity):

        totalgametocytes = 0

        if immunity:
            for i in range(GametocyteStages.Mature):
                totalgametocytes += self.m_malegametocytes[i] + self.m_femalegametocytes[i]

            totalIRBC = sum(self.m_IRBC_count)
            if totalIRBC > self.m_max_parasites:
                self.m_max_parasites = totalIRBC

            if (totalIRBC*self.m_inv_microliters_blood) > MINIMUM_IRBC_COUNT:
                self.m_start_measuring = True

            if self.m_start_measuring:
                self.m_temp_duration += dt
                if (totalIRBC*self.m_inv_microliters_blood) > MINIMUM_IRBC_COUNT:
                    self.m_measured_duration += self.m_temp_duration
                    self.m_temp_duration = 0

            if (totalIRBC+self.m_hepatocytes+totalgametocytes) < 1:
                self.StateChange = InfectionStateChange.Cleared

        return

    def SetParameters(self, strain=None, incubation_period_override=-1):

        self.CreateInfectionStrain(strain)

        if incubation_period_override != -1:
            self.incubation_timer = incubation_period_override
        else:
            self.incubation_timer = Config.incubation_distribution.CalculateDuration()

        self.infectious_timer = Config.infectious_distribution.CalculateDuration()
        self.total_duration = self.incubation_timer + self.infectious_timer
        self.infectiousness = 0
        self.StateChange = InfectionStateChange.NONE

        if self.incubation_timer <= 0:
            self.infectiousness = Config.base_infectivity

        return

    def CreateInfectionStrain(self, strain):

        if self.infection_strain is None:
            self.infection_strain = identity.StrainIdentity()

        if strain is not None:
            self.infection_strain.antigen_id = strain.antigen_id
            self.infection_strain.sub_strain = strain.sub_strain

        return
