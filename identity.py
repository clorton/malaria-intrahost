#!/usr/bin/python


class StrainIdentity(object):
    def __init__(self, antigen_id=0, sub_strain=0):

        self.antigen_id = antigen_id
        self.sub_strain = sub_strain

        return

    def ResolveInfectingStrain(self):

        return StrainIdentity(self.antigen_id, self.sub_strain)
