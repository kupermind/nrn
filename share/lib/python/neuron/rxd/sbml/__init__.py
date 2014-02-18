from neuron.rxd.rxdException import RxDException

try:
    import libsbml
except:
    raise RxDException('could not import libsbml; in Ubuntu use: "sudo apt-get install libsbml-python"')
try:
    import sympy
except:
    raise RxDException('could not import sympy; in Linux use: "sudo easy_install sympy"')

# TODO: relative imports don't work in Python 3 (I think)
from sbmlToRxd import load_sbml
import sbml_lists
