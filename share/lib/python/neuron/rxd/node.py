import neuron
from neuron import h, nrn
import rxd
import numpy
import weakref
import region
import rxdsection
from rxdException import RxDException

# data storage
_volumes = numpy.array([])
_surface_area = numpy.array([])
_diffs = numpy.array([])
_states = numpy.array([])

# node data types
_concentration_node = 0
_molecule_node = 1

def _get_data():
    return (_volumes, _surface_area, _diffs)

def _get_states():
    return _states

def _allocate(num):
    """allocate storage for num more nodes, return the starting index
    
    Note: no guarantee is made of preserving previous _ref
    """
    start_index = len(_volumes)
    total = start_index + num
    _volumes.resize(total, refcheck=False)
    _surface_area.resize(total, refcheck=False)
    _diffs.resize(total, refcheck=False)
    _states.resize(total, refcheck=False)
    return start_index

_numpy_element_ref = neuron.numpy_element_ref

class Node(object):
    def satisfies(self, condition):
        """Tests if a Node satisfies a given condition.

        If a nrn.Section object or RxDSection is provided, returns True if the Node lies in the section; else False.
        If a Region object is provided, returns True if the Node lies in the Region; else False.
        If a number between 0 and 1 is provided, returns True if the normalized position lies within the Node; else False.
        """
        if isinstance(condition, nrn.Section) or isinstance(condition, rxdsection.RxDSection):
            return self._in_sec(condition)
        elif isinstance(condition, region.Region):
            return self.region == condition
        try:
            dx = 1. / self.sec.nseg / 2.
            if 0 < condition <= 1:
                return -dx < self.x - condition <= dx
            elif condition == 0:
                # nodes at dx, 3dx, 5dx, 7dx, etc... so this allows for roundoff errors
                return self.x < 2. * dx
                
        except:
            raise RxDException('unrecognized node condition: %r' % condition)

    @property
    def _ref_concentration(self):
        """Returns a HOC reference to the Node's concentration
        
        (The node must store concentration data. Use _ref_molecules for nodes
        storing molecule counts.)
        """
        # this points to rxd array only, will not change legacy concentration
        if self._data_type == _concentration_node:        
            return self._ref_value
        else:
            raise RxDException('_ref_concentration only available for concentration nodes')

    @property
    def _ref_molecules(self):
        """Returns a HOC reference to the Node's concentration
        
        (The node must store concentration data. Use _ref_molecules for nodes
        storing molecule counts.)
        """
        # this points to rxd array only, will not change legacy concentration
        if self._data_type == _molecule_node:        
            return self._ref_value
        else:
            raise RxDException('_ref_molecules only available for molecule count nodes')
    
    @property
    def d(self):
        """Gets the diffusion rate within the compartment."""
        return _diffs[self._index]
    @d.setter
    def d(self, value):
        """Sets the diffusion rate within the compartment."""
        # TODO: make invalidation work so don't need to redo the setup each time
        #rxd._invalidate_matrices()
        _diffs[self._index] = value
        rxd._setup_matrices()    

    @property
    def concentration(self):
        """Gets the concentration at the Node."""
        # TODO: don't use an if statement here... put the if statement at node
        #       construction and change the actual function that is pointed to
        if self._data_type == _concentration_node:
            return _states[self._index]
        else:
            # TODO: make this return a concentration instead of an error
            raise RxDException('concentration property not yet supported for non-concentration nodes')
    
    @concentration.setter
    def concentration(self, value):
        """Sets the concentration at the Node"""
        # TODO: don't use an if statement here... put the if statement at node
        #       construction and change the actual function that is pointed to
        if self._data_type == _concentration_node:
            _states[self._index] = value
        else:
            # TODO: make this set a concentration instead of raise an error
            raise RxDException('concentration property not yet supported for non-concentration nodes')

    @property
    def molecules(self):
        """Gets the molecule count at the Node."""
        # TODO: don't use an if statement here... put the if statement at node
        #       construction and change the actual function that is pointed to
        if self._data_type == _molecule_node:
            return _states[self._index]
        else:
            # TODO: make this return a molecule count instead of an error
            raise RxDException('molecules property not yet supported for non-concentration nodes')
    
    @molecules.setter
    def molecules(self, value):
        """Sets the molecule count at the Node"""
        # TODO: don't use an if statement here... put the if statement at node
        #       construction and change the actual function that is pointed to
        if self._data_type == _molecule_node:
            _states[self._index] = value
        else:
            # TODO: make this set a molecule count instead of raise an error
            raise RxDException('molecules property not yet supported for non-concentration nodes')
        

        
    @property
    def value(self):
        """Gets the value associated with this Node."""
        return _states[self._index]
    
    @value.setter
    def value(self, v):
        """Sets the value associated with this Node.
        
        For Species nodes belonging to a deterministic simulation, this is a concentration.
        For Species nodes belonging to a stochastic simulation, this is the molecule count.
        """
        _states[self._index] = v

    @property
    def _ref_value(self):
        """Returns a HOC reference to the Node's value"""
        return _numpy_element_ref(_states, self._index)
    
_h_n3d = h.n3d
_h_x3d = h.x3d
_h_y3d = h.y3d
_h_z3d = h.z3d
_h_arc3d = h.arc3d

class Node1D(Node):
    def __init__(self, sec, i, location, data_type=_concentration_node):
        """n = Node1D(sec, i, location)
        Description:
        
        Constructs a Node object. These encapsulate multiple properties of a given
        reaction-diffusion compartment -- volume, surface area, concentration, etc... --
        but this association is provided only as a convenience to the user; all the data
        exists independently of the Node object.

        These objects are typically constructed by the reaction-diffusion code and
        provided to the user.

        Parameters:

        sec  -- the RxDSection containing the Node

        i -- the offset into the RxDSection's data

        location -- the location of the compartment within the section. For @aSection1D objects, this is the normalized position 0 <= x <= 1
        """
        self._sec = sec
        self._location = location
        self._index = i + sec._offset
        self._loc3d = None
        self._data_type = data_type


    def _update_loc3d(self):
        sec = self._sec
        length = sec.L
        normalized_arc3d = [_h_arc3d(i, sec=sec._sec) / length for i in xrange(int(_h_n3d(sec=sec._sec)))]
        x3d = [_h_x3d(i, sec=sec._sec) for i in xrange(int(_h_n3d(sec=sec._sec)))]
        y3d = [_h_y3d(i, sec=sec._sec) for i in xrange(int(_h_n3d(sec=sec._sec)))]
        z3d = [_h_z3d(i, sec=sec._sec) for i in xrange(int(_h_n3d(sec=sec._sec)))]
        loc1d = self._location
        self._loc3d = (numpy.interp(loc1d, normalized_arc3d, x3d),
                       numpy.interp(loc1d, normalized_arc3d, y3d),
                       numpy.interp(loc1d, normalized_arc3d, z3d))
        

    @property
    def x3d(self):
        """x coordinate"""
        if self._loc3d is None:
            self._update_loc3d()
        return self._loc3d[0]

    @property
    def y3d(self):
        """y coordinate"""
        if self._loc3d is None:
            self._update_loc3d()
        return self._loc3d[1]
    
    @property
    def z3d(self):
        """z coordinate"""
        if self._loc3d is None:
            self._update_loc3d()
        return self._loc3d[2]

    @property
    def volume(self):
        """The volume of the compartment in cubic microns.
        
        Read only."""
        rxd._update_node_data()
        return _volumes[self._index]
    
    @property
    def segment(self):
        return self._sec._sec(self.x)
    
    @property
    def surface_area(self):
        """The surface area of the compartment in square microns.
        
        This is the area (if any) of the compartment that lies on the plasma membrane
        and therefore is the area used to determine the contribution of currents (e.g. ina) from
        mod files or kschan to the compartment's concentration.
        
        Read only.
        """
        
        rxd._update_node_data()
        return _surface_area[self._index]
    
        
    @property
    def x(self):
        """The normalized position of the center of the compartment.
        
        Read only."""
        # TODO: will probably want to change this to be more generic for higher dimensions
        return self._location
    
    @property
    def region(self):
        """The region containing the compartment."""
        return self._sec._region
    
    @property
    def sec(self):
        """The RxDSection containing the compartment."""
        return self._sec

    def _in_sec(self, sec):
        return sec == self.sec or sec == self.sec._sec

    
    @property
    def species(self):
        """The Species whose concentration is recorded at this Node."""
        return self._sec._species()    


class Node3D(Node):
    def __init__(self, index, i, j, k, r, seg, speciesref, data_type=_concentration_node):
        """
            Parameters
            ----------
            
            index : int
                the offset into the global rxd data
            i : int
                the x coordinate in the region's matrix
            j : int
                the y coordinate in the region's matrix
            k : int
                the z coordinate in the region's matrix
            r : rxd.Region
                the region that contains this node
            seg : nrn.Segment
                the segment containing this node
        """
        self._index = index
        self._i = i
        self._j = j
        self._k = k
        self._r = r
        self._seg = seg
        self._speciesref = speciesref
        self._data_type = data_type
    
    @property
    def surface_area(self):
        """The surface area of the compartment in square microns.
        
        This is the area (if any) of the compartment that lies on the plasma membrane
        and therefore is the area used to determine the contribution of currents (e.g. ina) from
        mod files or kschan to the compartment's concentration.
        
        Read only.
        """
        # TODO: should I have the commented out line?
        #rxd._update_node_data()
        return _surface_area[self._index]
    
    @property
    def x3d(self):
        # TODO: need to modify this to work with 1d
        return self._r._mesh.xs[self._i]
    @property
    def y3d(self):
        # TODO: need to modify this to work with 1d
        return self._r._mesh.ys[self._j]
    @property
    def z3d(self):
        # TODO: need to modify this to work with 1d
        return self._r._mesh.zs[self._k]
    
    @property
    def x(self):
        raise RxDException('need to reimplement x for 3d nodes')
    
    @property
    def segment(self):
        return self._seg
    
    def _in_sec(self, sec):
        return sec == self.sec
    
    @property
    def sec(self):
        if self._seg is None:
            return None
        return self._seg.sec
    
    @property
    def volume(self):
        return _volumes[self._index]

    @property
    def region(self):
        """The region containing the compartment."""
        return self._r

    @property
    def species(self):
        """The Species whose concentration is recorded at this Node."""
        return self._speciesref()
