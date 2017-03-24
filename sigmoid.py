#!/usr/bin/python


def basic_sigmoid(threshold = 100.0, variable = 0.0):
    return (variable/(threshold+variable)) if variable > 0 else 0.0
