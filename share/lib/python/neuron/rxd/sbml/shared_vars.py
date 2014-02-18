''' shared_vars.py contains variable shared by more than one file in the package '''

#globals
species_amount_conc = {}
scaling = {}

#sbml main objects
compartments = {}
species = {}
species_keys = [] #replaced by speciesArray
parameters = {}
parameters_keys = []
derivs = {}
rate_eqs = {}
species_derivs = {}
assignment_rules = {}
assignment_keys = []
reactions = {}


# arrays
regions = []
speciesArray = {}

#rxd lists
rxd_compartments = {}
rxd_species = {}
rxd_parameters = {}
rxd_rates = {}
rxd_derivs = {}

#events 
eventsTrigs = {}
eventsDels = {}
eventsAssigs = {}
