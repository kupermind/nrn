#!/usr/bin/env python
# coding=utf-8
''' sbml_numpy.py contains math functions to be used in sbml to rxd (neuron) converter '''

import sys, os, libsbml
from math import *
from operator import *
from sympy import symbols
from sympy.physics.units import *
from decimal import *
import xml.dom.minidom

from neuron import h, rxd
from numpy import *
from matplotlib import pyplot
import shlex


def power_safe(a, b):
	try:
		return float(a) ** float(b)
	except ValueError:
		return -(abs(a) ** float(b))

def arccosh(x):
	a = math.sqrt(x**2 - 1)
	if x <= 0:
		return math.log(x - a)
	else:
		return math.log(x + a)

def arctanh(x):
	if math.fabs(x) > 1:
		raise ValueError("Area tangens hyperbolicus is defined in the interval [-1, 1] only")
	return math.log(math.sqrt( (1+x)/(1-x) ))

def arccoth(x):
	if math.fabs(x) <= 1:
		raise ValueError("Area cotangens hyperbolicus is not defined in the interval [-1, 1]")
	return math.log(math.sqrt((x+1)/(x-1)))

def arcsech(x):
	a = math.sqrt(x**2 - 1)
	if x <= 0:
		return (1 / math.log(x - a))
	else:
		return (1 / math.log(x + a))

sin = math.sin
cos = math.cos
cot = lambda x: math.cos(x) / math.sin(x)
sec = lambda x: 1 / math.cos(x)
csc = lambda x: 1 / math.sin(x)
coth = lambda x: math.cosh(x) / math.sinh(x)
sech = lambda x: 2 / (math.exp(x) - math.exp(-x))
csch = lambda x: 2 / (math.exp(x) - math.exp(-x))
arcsin = math.asin
arccos = math.acos
arctan = math.atan
arccot = lambda x: math.acos(x / math.sqrt(1 + x**2))
arcsec = lambda x: 1 / math.acos(x)
arccsc = lambda x: 1 / math.asin(x)
arcsinh = lambda x: math.log(x + math.sqrt(x**2  + 1))
arccosh = arccosh
arctanh = arctanh
arccoth = arccoth
arcsech = arcsech
arccsch = lambda x: 1/math.log(x + math.sqrt(x**2 + 1))
ln = math.log
power = power_safe
root = lambda a,b: math.pow(a, 1/b)
exp = lambda x: math.exp(x)
factorial = lambda x: reduce(operator.mul, xrange(1, x + 1))
log = math.log

qeg = ge
leq = le
neq = ne
sbml_and = lambda a, b: bool(a) and bool(b)
sbml_not = lambda a: not bool(a)
sbml_or = lambda a, b: bool(a) and bool(b) 
xor = lambda a, b: ((not a and b) or ( a and not b))


# TYPES TO BE HANDLED
#   AST_NAME_TIME = _libsbml.AST_NAME_TIME
#	AST_LAMBDA = _libsbml.AST_LAMBDA
#	AST_FUNCTION = _libsbml.AST_FUNCTION
#	AST_FUNCTION_DELAY = _libsbml.AST_FUNCTION_DELAY
#	AST_FUNCTION_PIECEWISE = _libsbml.AST_FUNCTION_PIECEWISE


# Reconcile american and british english
molar = mol/liter
litre = liter
metre = meter
item, items = symbols("item items")
item = mole/avogadro_number
items = mole/avogadro_number

sympy_units = {
	"volume": meter**3,
	"area": meter**2,
	"length": meter
}

sbml_units = {
	"time": second, 
    "substance": mole,		
    "volume": liter,	
    "area": meter**2,	
    "length": meter,		
}

rxd_units = {
	"time": milli*second,
	"concentration": milli*(mole/liter),
	"volume": micron**3,
	"area": micron**2,
	"length": micron,
	"mV": milli*volt,
    "mA": milli*ampere
}

abortSimulationList = []
def recombineString(split_formula, oldTerm, newTerm, brackets):
	new_formula = ''

	for token in split_formula:
		if token == "piecewise":
			errorMsg = "a formula in this model uses piecewise function which is not currently supported"
			if errorMsg not in abortSimulationList:
				abortSimulationList.append(errorMsg)

		if token == oldTerm:
			if brackets == True:
				new_formula += "(" +newTerm+ ")"
			else:
				new_formula += "" +newTerm+ ""
		else:
			new_formula += token

	return new_formula


def replaceInFormula(formula, oldTerm, newTerm, brackets=True):
	split_formula = list(shlex.shlex(formula))
	return recombineString(split_formula, oldTerm, newTerm, brackets)


def replaceNegValues(formula, oldTerm, newTerm, brackets=False):
	split_formula = list(shlex.shlex(formula))

	if split_formula[0] == "-":
		split_formula[0] = "-1.0*"

	return recombineString(split_formula, oldTerm, newTerm, brackets)



