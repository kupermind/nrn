#!/usr/bin/env python


# 	filename: sbmlToRxd.py
# 	brief: it extracts data from SBML files and uses it to run reaction-diffusion simulations using neuron.rxd
# 	author: Victor Mutai, muvic08@gmail.com
# 	references: 
#		SBML:   http://sbml.org/Software/libSBML/docs/java-api/org/sbml/libsbml/Model.html
#				http://www.ebi.ac.uk/biomodels/models-main/publ/BIOMD0000000340/BIOMD0000000340.pdf
# 		RXD: 	http://www.neuron.yale.edu/neuron/static/new_doc/modelspec/programmatic/rxd.html
#				http://www.neuron.yale.edu/neuron/static/docs/rxd/index.html
#		CREDITS: Padraig Gleeson for SBML2NEURON.py which informed my earlier decision


# import local modules
from sbml_math import *
from sbml_lists import *
import imp


def printToFile(outfile, rxdfile):
	rxdfile.write("#!/usr/bin/env python\n\nfrom neuron import h, rxd \nfrom matplotlib import pyplot \n")
	rxdfile.write("import sys\nsys.path.append('{}')\nfrom sbml_math import *\nfrom operator import *\n\n".format(os.getcwd()))
	rxdfile.write('h.load_file("stdrun.hoc")\nrxd.options.use_reaction_contribution_to_jacobian = False\n\n')

	rxdfile.write("soma = h.Section() \ndend = h.Section() \ntime = rxd.Parameter(rxd.Region(h.allsec()), initial=0.0) \nt = time \n")

	rxdfile.write("\n#PARAMETERS:\n")
	rxdfile.write("parameters_keys = {}\n\n".format(parameters_keys))
	for key in rxd_parameters.keys():
		rxdfile.write("{} = {}\n".format(key, rxd_parameters[key]))

	rxdfile.write("\n#RULES:\n")
	rxdfile.write("assignment_keys = {}\n".format(assignment_keys))
	rxdfile.write("assignment_rules = {}\n".format(assignment_rules))

	rxdfile.write("\n#COMPARTMENTS:\n")
	for key in rxd_compartments.keys():
		rxdfile.write("{} = {}\n".format(key, rxd_compartments[key]))

	rxdfile.write("\n#SPECIES:\n")
	rxdfile.write("speciesArray = {}\n\n".format(rxd_species.keys()))
	for key in rxd_species.keys():
		rxdfile.write("{} = {}\n".format(key, rxd_species[key]))

	rxdfile.write("\n#DERIVATIVES:\n")
	for key in rxd_derivs.keys():
		rxdfile.write("{} = {}\n".format(key, rxd_derivs[key]))

	rxdfile.write("\n#RATES:\n")
	for key in rxd_rates.keys():
		rxdfile.write("{} = {}\n".format(key, rxd_rates[key]))

	rxdfile.write("\n#EVENTS:\n")
	rxdfile.write("eventsTrigs = {}\n".format(eventsTrigs))
	rxdfile.write("eventsAssigs = {}\n".format(eventsAssigs))
	rxdfile.write("eventsDels = {}\n".format(eventsDels))

	rxdfile.write("\n#ERRORS:\nsimulationErrors = {}\n\n".format(abortSimulationList))


	rxdfile.write("")


	outfile.write("\nUNITS: \n\n")
	for unit in sbml_units:
		outfile.write("{} = {}\n".format(unit, sbml_units[unit]))

	outfile.write("\n\nCOMPARTMENTS: \n\n")
	for c in compartments:
		outfile.write("{} = {}\n".format(c, compartments[c]))

	outfile.write("\n\nSPECIES: \n\n")
	for s in species:
		outfile.write("{} = {}\n".format(s, species[s]))

	outfile.write("\n\nPARAMETERS: \n\n")
	for p in parameters:
		outfile.write("{} = {}\n".format(p, parameters[p]))	

	outfile.write("\n\nREACTIONS: \n\n")
	for r in reactions:
		outfile.write("{} = {}\n".format(r, reactions[r]))

	outfile.write("\n\nDERIVS")
	for d in species_derivs:
		outfile.write("{} = {}\n".format(d, species_derivs[d]))



def rxdReactions():
	for spes in species:
		ref = spes+"_rate"
		der = spes+"_deriv"

		if spes in derivs: 
			if speciesArray[spes]['constant'] is False and speciesArray[spes]['boundaryCondition'] is False:
				rxd_derivs[der] = derivs[spes]
				rxd_rates[ref] = "rxd.Rate("+spes+", eval('"+der+"'))"

		if spes in rate_eqs:
			if speciesArray[spes]['constant'] is False and speciesArray[spes]['boundaryCondition'] is True:
				rxd_derivs[der] = rate_eqs[spes]
				rxd_rates[ref] = "rxd.Rate("+spes+", eval('"+der+"'))"



def getLists(theList, theName, printFunc, outfile):
	mod = sys.modules['__main__']

	printFuction = getattr(mod, printFunc) 
	outfile.write("#	{} {}\n".format(len(theList), theName))
	for n, l in enumerate(theList):
		printFuction(n+1, l)




def getListsFromSBML(model, infile, outfile):
	outfile.write("#SBML file ({}) has: \n".format(infile))

	sbmlModelAttr = [
		"UnitDefinitions", 
		"Compartments", 
		"Species", 
		"Parameters", 
		"Rules", 
		"Constraints", 
		"Reactions",
		"Events"
	]

	def getLists(theList, theName, printFunc, outfile):
		mod = sys.modules['__main__']

		printFuction = getattr(mod, printFunc) 
		outfile.write("#	{} {}\n".format(len(theList), theName))
		for l in theList:
			printFuction(l)


	for constant in sbmlModelAttr:
		theList = "getListOf"+constant
		func = getattr(model, theList)
		if constant is not "Species":
			printFunc = "get"+constant[:-1]
		else:
			printFunc = "get"+constant

		if func():
			getLists(func(), constant, printFunc, outfile)


def runSimulations(raw_rxdfile_dir, file_number):
	raw_filename = "BIOMD0000000"+file_number+".py"
	generatedFileStatement = "You can run (or edit) simulations on the generated file ({}). Full file path below: \n\t{} ".format(raw_filename, raw_rxdfile_dir+"/"+raw_filename)
	print generatedFileStatement+"\n"

	#TODO: run the newly generated file directly from here

	#runSim = raw_input("Do you want to run {} right away? [y] [n] ".format(raw_filename+".py"))
	#if runSim is not "n" or "N":
	#	sys.path.append(raw_rxdfile_dir)
	#	module = __import__(raw_filename)
	#else:
	#	print "You can  " '''

def main(args):

	# define global variables
	global assignment_rules_counter
	assignment_rules_counter = True 

	# directories and files
	modelsDir = os.getcwd()
	biomodelsDir = "/biomodels-2013/curated"
	commonFileStr = "/BIOMD0000000"
	
	# create a directories to store output data
	if 'biomodels-output' in os.listdir(os.curdir): pass
	else: os.mkdir("biomodels-output")

	if 'biomodels-rxd' in os.listdir(os.curdir): pass
	else: os.mkdir("biomodels-rxd")

	file_number = raw_input("\nEnter the last 3 digits of your biomodel (BIOMD0000000***) filename [001 - 462]: ")

	while len(file_number) != 3 or int(file_number) > 462 or int(file_number) < 0:
		file_number = raw_input("Please re-enter the last 3 digits of your biomodel filename [001 - 462]: ")

	infile = os.path.abspath(modelsDir+biomodelsDir+commonFileStr+file_number+".xml")
	outfile = os.path.abspath(modelsDir+"/biomodels-output/"+commonFileStr+file_number+".txt")
	rxdfile = os.path.abspath(modelsDir+"/biomodels-rxd/"+commonFileStr+file_number+".py")
	raw_rxdfile_dir = os.path.abspath(modelsDir+"/biomodels-rxd/")
	
	print "Biomodel file path: ", infile, "\n"

	outfile = file(outfile, mode='w')
	rxdfile = file(rxdfile, mode='w')

	domFile = file(infile, mode="r")

	dom = xml.dom.minidom.parseString(domFile.read())

	species = dom.getElementsByTagName("species")	

	for specis in species:
		if specis.hasAttribute("initialAmount"):
			species_amount_conc["conc_"+specis.getAttribute("id")] = "initialAmount"
		elif specis.hasAttribute("initialConcentration"):
			species_amount_conc["conc_"+specis.getAttribute("id")] = "initialConcentration"
		else:
			species_amount_conc["conc_"+specis.getAttribute("id")] = "unknown"

	domFile.close()

   	# read sbml from infile
   	reader = libsbml.SBMLReader()
	sbmldoc = reader.readSBML(infile)

	# expand function defininitions and initial assignments
	libsbml.SBMLDocument.expandFunctionDefinitions(sbmldoc)
	libsbml.SBMLDocument.expandInitialAssignments(sbmldoc)

	# get sbml model and extract the data
	model = sbmldoc.getModel()
	getListsFromSBML(model, infile, outfile)
	
	rxdReactions() # create a list of derivatives 

	if len(speciesArray) < 1:
		abortSimulationList.append("this model has 0 species")
		
	if len(rxd_compartments) > 1:
		abortSimulationList.append("this model has more than one compartment ("+str(len(rxd_compartments))+")")
		

	runSimulations(raw_rxdfile_dir, file_number)


	printToFile(outfile, rxdfile)

	with open("rxd_simulation.py", "r") as simFile:
		data = simFile.read()
		rxdfile.writelines(data)




if __name__ == '__main__':
    main(sys.argv)
