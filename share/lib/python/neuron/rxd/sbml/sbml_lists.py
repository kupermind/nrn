''' sbml_lists.py extracts data from sbml. This file includes the following sbml properties in the simulations:
		"UnitDefinitions", 
		"Compartments", 
		"Species", 
		"Parameters", 
		"Rules", 
		"Reactions",
		"Events"
'''


from sbml_math import *
from shared_vars import *

def getUnitDefinition(listItem):
	theId = listItem.getId()
	for unit in listItem.getListOfUnits():
		kind = unit.getKind()
		kindToString = libsbml.UnitKind_toString(kind)
		multiplier = unit.getMultiplier()
		scale = unit.getScale()
		exponent = unit.getExponent()
		offset = unit.getOffset()		

		if kindToString != "dimensionless":
			newUnit = (multiplier*(10**scale)*eval(kindToString))**exponent + offset
			sbml_units[theId] = newUnit




def	getCompartment(listItem):
	theId = listItem.getId()
	regions.append(theId)

	spatial_dimension = listItem.getSpatialDimensions()
	global space

	if spatial_dimension == 1: space = "length"
	elif spatial_dimension == 2: space = "area"
	elif spatial_dimension == 3: space = "volume"

	spatial_value = sbml_units[space]/sympy_units[space] 

	compartments[theId] = listItem.getSize()
	globals()[theId] = compartments[theId]

	rxd_compartments[theId] = "rxd.Region(h.allsec())"




def getSpecies(listItem):
	theId = listItem.getId()
	species_keys.append(theId)

	if species_amount_conc["conc_"+listItem.getId()] == "initialConcentration" or species_amount_conc["conc_"+listItem.getId()] == "unknown":
		initial = listItem.getInitialConcentration()
		species_units = sbml_units["substance"] / sbml_units["volume"]
	elif species_amount_conc["conc_"+listItem.getId()] == "initialAmount": 
		# Convert amount into concentration by dividing it with the volume of the compartment 
		initial = listItem.getInitialAmount()
		species_units = sbml_units["substance"] / (compartments[listItem.getCompartment()] * sbml_units["volume"])
	else:
		print "I'VE NO IDEA WHAT IS HAPPENING! CHECK THE SPECIES"

	species[theId] = float(initial)
	speciesObject = {"constant": listItem.getConstant(), "boundaryCondition": listItem.getBoundaryCondition()}
	speciesArray[theId] = speciesObject

	scaling_value = (species_units)/(milli*molar)
	scaling[listItem.getId()] = scaling_value

	initialAmount_in_molar = initial*scaling_value  
	initialAmount_in_str = "%g" % initialAmount_in_molar
	rxd_species[listItem.getId()] = "rxd.Species("+listItem.getCompartment()+", d=1, initial=float("+str(initialAmount_in_str)+"))"




def getParameter(listItem):
	parameters[listItem.getId()] = listItem.getValue()
	parameters_keys.append(listItem.getId())
	rxd_parameters[listItem.getId()] = "rxd.Parameter(rxd.Region(h.allsec()), initial="+str(listItem.getValue())+")" 




def getRule(listItem):
	if listItem.getType() == libsbml.RULE_TYPE_RATE:
		variable = listItem.getVariable()
		formula = listItem.getFormula()

		if variable in rate_eqs.keys(): 
			rate_eqs[variable] += "+("+formula+")"
		else: 
			rate_eqs[variable] = "+("+formula+")"

		if variable in species_derivs.keys():
			species_derivs[variable] += "+("+formula+")"
		else:
			species_derivs[variable] = "+("+formula+")"



	if listItem.getType() == libsbml.RULE_TYPE_SCALAR:
		formula = listItem.getFormula()
		assignment_rules[listItem.getVariable()] = formula
		assignment_keys.append(listItem.getVariable())




def getConstraints(listItem):
	pass

def replaceSpecies(formula):
		for spec in species_keys:
			#newValue =  ("((1.0/(%f)" % scaling[spec])+") * eval('"+str(spec)+"').nodes[0].concentration)"
			scaler = "1.0/(%g)" % scaling[spec]
			newValue =  str(spec)+"*("+scaler+")" #  ("(1.0/(%g))" % scaling[spec])+")"+str(spec)+")"
			formula = replaceInFormula(formula, spec, newValue, False)

		return formula

def getReaction(listItem):
	reactants = {}
	reactants_keys = []
	reactScheme = ""
	products = {}
	products_keys = []
	prodScheme = ""
	reactionParams = {}
	formula = ""
	oFormula = ""

	if listItem.isSetKineticLaw():
		kinecticLaw = listItem.getKineticLaw()
		if kinecticLaw.isSetMath():
			formula = kinecticLaw.getFormula()
			oFormula = formula

		if kinecticLaw.getListOfParameters():
			for p in kinecticLaw.getListOfParameters():
				reactionParams[p.getId()] = p.getValue()



	if listItem.getListOfProducts():
		for product in listItem.getListOfProducts():
			if product.getStoichiometry() != 1:
				molecule = str(product.getStoichiometry()) + "*" + str(product.getSpecies())
			else:
				molecule = str(product.getSpecies())

			products[product.getSpecies()] = molecule
			#newValue =  ("((1.0/(%g)" % scaling[product.getSpecies()])+") * "+str(product.getSpecies())+")"
			#formula = replaceInFormula(formula, str(product.getSpecies()), newValue)

			if prodScheme != "":
				prodScheme = "{} + {}".format(prodScheme, molecule)
			else:
				prodScheme = "{}".format(molecule)

			reaction_scaler = scaling[product.getSpecies()]
			reactants_keys.append(product.getSpecies)

	if listItem.getListOfReactants():
		for reactant in listItem.getListOfReactants():
			if reactant.getStoichiometry() != 1:
				molecule = str(reactant.getStoichiometry()) + "*" + str(reactant.getSpecies())
			else:
				molecule = str(reactant.getSpecies())

			reactants[reactant.getSpecies()] = molecule
			#newValue =  ("((1.0/(%g)" % scaling[reactant.getSpecies()])+") * "+str(reactant.getSpecies())+")"
			#formula = replaceInFormula(formula, str(reactant.getSpecies()), newValue)

			if reactScheme != "":
				reactScheme = "{} + {}".format(reactScheme, molecule)
			else:
				reactScheme = "{}".format(molecule)

			reaction_scaler = scaling[reactant.getSpecies()]
			products_keys.append(reactant.getSpecies())

	if listItem.getReversible():
		reactionArrow = "<>"
		rate1 = 1.0
		rate2 = 1.0
	else:
		reactionArrow = ">"
		rate1 = 1.0
		rate2 = None

	if len(reactScheme) == 0:
		reactScheme = '0 * dummy_species'

	if len(prodScheme) == 0:
		prodScheme = '0 * dummy_species'


	scheme = "{} {} {}".format(reactScheme, reactionArrow, prodScheme)

	for p in reactionParams:
		formula = replaceInFormula(formula, str(p), str(reactionParams[p]))

	for comp in compartments:
		newValue = "%g" % 1.0
		formula = replaceInFormula(formula, str(comp), str(newValue))

	timeRatio = sbml_units["time"] / rxd_units["time"]
	formula = replaceNegValues(formula, "", "")
	
	formula = replaceSpecies(formula)

	finalFormula = "(%s)*(%g)/(%g)" %(formula, float(reaction_scaler), timeRatio)
	#finalFormula = "(%s)*(%g)" %(formula, float(reaction_scaler))

	reactions[listItem.getId()] = scheme + "\n\t"+finalFormula 

	globals()[listItem.getId()] = finalFormula

	for prod in products:
		if prod in derivs.keys():
			species_derivs[prod] += "+("+oFormula+")"
			derivs[prod] += "+1.0*("+finalFormula+")"
		else:
			species_derivs[prod] = "+("+oFormula+")"
			derivs[prod] = "+1.0*("+finalFormula+")"

	for react in reactants:
		if react in derivs.keys():
			species_derivs[react] += "-("+oFormula+")"
			derivs[react] += "-1.0*("+finalFormula+")"
		else:
			species_derivs[react] = "-("+oFormula+")"
			derivs[react] = "-1.0*("+finalFormula+")"


def replaceSpeciesParamsRegs(formula):
		for spec in species_keys:
			#newValue =  ("((1.0/(%f)" % scaling[spec])+") * eval('"+str(spec)+"').nodes[0].concentration)"
			newValue =  ("((1.0/(%g)" % scaling[spec])+") * "+str(spec)+".nodes[0].concentration)"
			formula = replaceInFormula(formula, spec, newValue, False)

		for param in parameters_keys:
			#newValue = "eval('"+str(param)+"').nodes[0].concentration"
			newValue = str(param)+".nodes[0].concentration"
			formula = replaceInFormula(formula, param, newValue, False)

		for region in regions:
			formula = replaceInFormula(formula, region, str(compartments[region]), False)

		return formula


def getEvent(listItem):
	assignmentsList = []

	def gt_time(A, B):
		timeRatio = float(sbml_units["time"] / rxd_units["time"])
		B *= timeRatio
		return B


	if listItem.isSetTrigger():
		formula = libsbml.formulaToString(listItem.getTrigger().getMath()) or libsbml.formulaToString(listItem.getDelay())
		#formula = replaceInFormula(formula, "gt", "gt_time", False)
		formula = replaceSpeciesParamsRegs(formula)

		sbml_logical_operators = ['and', 'or', 'not']

		for op in sbml_logical_operators:
			formula = replaceInFormula(formula, op, "sbml_"+str(op), False)

		eventsTrigs[listItem.getId()] = formula
		

	for n, ea in enumerate(listItem.getListOfEventAssignments()):
		if ea.isSetMath():
			eaObj = { "variable": ea.getVariable(), "formula": replaceSpeciesParamsRegs(libsbml.formulaToString(ea.getMath())) }
			assignmentsList.append(eaObj)

	eventsAssigs[listItem.getId()] = assignmentsList

	eventsDels[listItem.getId()] = "0"
	if listItem.isSetDelay():
		formula = libsbml.formulaToString(listItem.getDelay().getMath()) or libsbml.formulaToString(listItem.getDelay())
		eventsDels[listItem.getId()] = formula

