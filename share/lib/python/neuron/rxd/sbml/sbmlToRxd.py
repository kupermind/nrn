#!/usr/bin/env python


#     filename: sbmlToRxd.py
#     brief: it extracts data from SBML files and uses it to run reaction-diffusion simulations using neuron.rxd
#     author: Victor Mutai, muvic08@gmail.com
#     references: 
#        SBML:   http://sbml.org/Software/libSBML/docs/java-api/org/sbml/libsbml/Model.html
#                http://www.ebi.ac.uk/biomodels/models-main/publ/BIOMD0000000340/BIOMD0000000340.pdf
#         RXD:     http://www.neuron.yale.edu/neuron/static/new_doc/modelspec/programmatic/rxd.html
#                http://www.neuron.yale.edu/neuron/static/docs/rxd/index.html
#        CREDITS: Padraig Gleeson for SBML2NEURON.py which informed my earlier decision

import tempfile


# import local modules
from sbml_math import *
from sbml_lists import *

# TODO: is this really needed? why?
import imp

def printToFile(outfile, rxdfile):
    rxdfile.write("#!/usr/bin/env python\n\nfrom neuron import h, rxd\n")
    # TODO: probably don't need the path changing stuff
    rxdfile.write("import sys\nsys.path.append('{}')\nfrom neuron.rxd.sbml.sbml_math import *\nfrom operator import *\n\n".format(os.getcwd()))
    rxdfile.write('h.load_file("stdrun.hoc")\nrxd.options.use_reaction_contribution_to_jacobian = False\n\n')


    # TODO: what is the 't' or 'time' parameter for?
    #       I'm guessing this should provide a way of accessing h.t (it doesn't now)
    rxdfile.write("time = rxd.Parameter(rxd.Region(h.allsec()), initial=0.0) \nt = time \n")

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
    outfile.write("#    {} {}\n".format(len(theList), theName))
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
        # TODO: this is just a way of selecting functions from the module; do this more directly!
        mod = sys.modules['neuron.rxd.sbml.sbml_lists']

        printFuction = getattr(mod, printFunc) 
        outfile.write("#    {} {}\n".format(len(theList), theName))
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


rxd_file_boilerplate = """
if len(simulationErrors) > 0:
    print "Simulation aborted because:"
    for err in simulationErrors:
        print "\t", err
    from neuron.rxd.rxdException import RxDException
    RxDException('sbml import issues')

def simRule(rule_key, sim_time=0):
    formula = assignment_rules[rule_key]
    for key in parameters_keys:
        formula = replaceInFormula(formula, str(key), str(globals()[key].nodes[0].concentration))

    for key in speciesArray:
        formula = replaceInFormula(formula, str(key), str(globals()[key].nodes[0].concentration))

    globals()[rule_key].nodes[0].concentration = eval(formula)


def get_sbml_time(rxd_time):
    timeRatio = float(sbml_units["time"] / rxd_units["time"])
    sbml_time = rxd_time / timeRatio
    return sbml_time #rxd_time


def assignValues(key):
    for event in eventsAssigs[key]:
        node = eval(event["variable"]).nodes[0]
        node.value = float(eval(event["formula"]))


h.finitialize()
y_values = []
t_values = []
assignment_rules_counter = True

if len(eventsTrigs):
    from neuron.rxd.rxdException import RxDException
    RxDException('events currently not supported, but see original version from Victor')

"""


sbml_import_count = 0

def load_sbml(infile):
    global sbml_import_count
    # define global variables
    global assignment_rules_counter
    assignment_rules_counter = True 
    
    sbml_import_count += 1
    
    # TODO: this is bad; could overwrite existing file or be edited maliciously by another program; also ugly
    rxdfile_name = 'tmp_sbml_%d' % sbml_import_count

    # TODO: remove outfile... might keep for now though, because provides diagnostic info
    outfile = tempfile.TemporaryFile()
    rxdfile = open(rxdfile_name + '.py', 'w')

    with open(infile) as domFile:
        dom = xml.dom.minidom.parseString(domFile.read())

    species = dom.getElementsByTagName("species")    

    for specis in species:
        if specis.hasAttribute("initialAmount"):
            species_amount_conc["conc_"+specis.getAttribute("id")] = "initialAmount"
        elif specis.hasAttribute("initialConcentration"):
            species_amount_conc["conc_"+specis.getAttribute("id")] = "initialConcentration"
        else:
            species_amount_conc["conc_"+specis.getAttribute("id")] = "unknown"

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
        
    printToFile(outfile, rxdfile)
    rxdfile.write(rxd_file_boilerplate)
    rxdfile.close()
    
    exec('import ' + rxdfile_name)
    os.remove(rxdfile_name + '.py')
    
    return locals()[rxdfile_name]

