if len(simulationErrors) > 0:
	print "Simulation aborted because:"
	for err in simulationErrors:
		print "\t", err
	sys.exit("\n")


#RUN SIMULATION:

print "\nspeciesArray:", speciesArray

species_name = raw_input("\nInput one species in the speciesArray above: ")
while species_name not in speciesArray:
	species_name = raw_input("Invalid input: enter the the species name from the array: ")

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

lastState = {}
for key in eventsTrigs.keys():
	lastState[key] = False


def plot_it(color='k'):
	time = get_sbml_time(h.t)
	t = time

	counter = 0
	timesObj = {}
	for key in assignment_keys:
		simRule(key, h.t)
		'''if globals()["assignment_rules_counter"]:
			counter += 1
			simRule(key, h.t)
			#print "updating assignm value of:", key
			#if counter == len(assignment_keys):
			print "\n"
			globals()["assignment_rules_counter"] = False'''

	for key in eventsTrigs.keys():
		def eventFuction():
			assignValues(key)

		if eval(eventsTrigs[key]) == True and h.t < 1:
			lastState[key] = True

		if eval(eventsTrigs[key]) == True and lastState[key] == False:
			print "key", key
			print "time:", time
			if eval(eventsDels[key]) is 0:
				assignValues(key)
			else:
				h.CVode().event(time+eval(eventsDels[key]), eventFuction)
			lastState[key] = True

		elif eval(eventsTrigs[key]) == False:
			lastState[key] = False
		
	simObj = eval(species_name)  
	y = simObj.nodes[0].concentration
	x = simObj.nodes[0].x
	print species_name, h.t, y	

	y_values.append(y)
	t_values.append(h.t)

plot_it('r')

for i in xrange(1, 10000, 1):
	h.CVode().event(i*.1, plot_it)

h.dt = .025
h.CVode().active(1)
h.continuerun(200)
pyplot.plot(t_values, y_values, 'r')
pyplot.show()
