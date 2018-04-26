# Converts the SBLI 2d block grid into HDF5 format. For reading into opensbli
import numpy as np
import h5py
#from opensbli import *
#from opensbli.utilities.helperfunctions import set_hdf5_metadata

class SimulationBlock(object):  # BoundaryConditionTypes add this later
    """
    """

    def __init__(self, ndim, block_number=None):
        if block_number:
            self.__blocknumber = block_number
        else:
            self.__blocknumber = 0
        self.ndim = ndim
        #KernelCounter.__init__(self)
        #Grid.__init__(self)
        # RationalCounter.__init__(self)
        self.boundary_halos = [[set(), set()] for d in range(self.ndim)]
        self.block_datasets = {}
        self.constants = {}
        self.Rational_constants = {}
        self.block_stencils = {}
        self.InputOutput = []
        self.list_of_equation_classes = []
        return

    @property
    def blockname(self):
        """
        """
        return 'opensbliblock%02d' % self.blocknumber

    @property
    def blocknumber(self):
        return self.__blocknumber

    @blocknumber.setter
    def blocknumber(self, number):
        self.__blocknumber = number
        return


def set_hdf5_metadata(dset, halos, npoints, block):
    """ Function to set hdf5 metadata required by OPS to a dataset. """
    d_m = [halos[0]]*block.ndim
    d_p = [halos[1]]*block.ndim

    dset.attrs.create("d_p", d_p, dtype="int32")
    dset.attrs.create("d_m", d_m, dtype="int32")
    dset.attrs.create("dim", [1], dtype="int32")
    dset.attrs.create("ops_type", u"ops_dat",dtype="S7")
    dset.attrs.create("block_index", [block.blocknumber], dtype="int32")
    dset.attrs.create("base", [0 for i in range(block.ndim)], dtype="int32")
    dset.attrs.create("type", u"double",dtype="S15")
    dset.attrs.create("block", u"%s" % block.blockname,dtype="S25")
    dset.attrs.create("size", npoints, dtype="int32")
    return

def apply_group_attributes(group, block):
    group.attrs.create("dims", [block.ndim], dtype="int32")
    group.attrs.create("ops_type", u"ops_block",dtype="S8")
    group.attrs.create("index", [block.blocknumber], dtype="int32")
    return

def output_hdf5(array, array_name, halos, npoints, block):
    """ Creates an HDF5 file for reading in data to a simulation, 
    sets the metadata required by the OPS library. """
    if not isinstance(array, list):
        array = [array]
    if not isinstance(array_name, list):
        array_name = [array_name]
    assert len(array) == len(array_name)
    with h5py.File('data.h5', 'w') as hf:
        # Create a group
        if (isinstance(block, MultiBlock)):
            all_blocks = block.blocks
        else:
            all_blocks = [block]
        for b in all_blocks:
            g1 = hf.create_group(b.blockname)
            # Loop over all the dataset inputs and write to the hdf5 file
            for ar, name in zip(array, array_name):
                g1.attrs.create("dims", [b.ndim], dtype="int32")
                g1.attrs.create("ops_type", u"ops_block",dtype="S9")
                g1.attrs.create("index", [b.blocknumber], dtype="int32")
                block_dset_name = b.location_dataset(name).base
                dset = g1.create_dataset('%s' % (block_dset_name), data=ar)
                set_hdf5_metadata(dset, halos, npoints, b)
    return

block_files = ["Bl1.dat", "Bl2.dat","Bl3.dat"]
#block_files = ["Bl2.dat"]
#h5f = h5py.File('bl2.h5', 'w')
nhalo = 5
stepx = 1
stepy = 1
fname = "data1.h5"
h5f = h5py.File(fname, 'w')
for no, b in enumerate(block_files):
    print b
    f=open(b)
    nx,ny = map(int, f.readlines()[0].split())
    print nx, ny
    f.close()
    x,y,z = np.loadtxt(b, skiprows =1,unpack=True)
    x = x.reshape(nx, ny)
    y = z.reshape(nx, ny)
    # Take away the corner point from block 0 and 2
    #if no != 1:
        #x = x[1:,:]
        #y = y[1:,:]
    
    if no == 1:
        x = x[1:-1,:]
        y = y[1:-1,:]
    # For the first block rotate the data so that left corner becomes the starting 
    #print x[0,0], x[-1,0]
    #if no == 0:
        #x = x[::-1,:]
        #y = y[::-1,:]
    #print x[0,0], x[-1,0]
    shape = x.shape
    print shape
    new_shape = tuple([s+2*nhalo for s in reversed(shape)])
    print new_shape
    newx = np.zeros(new_shape)
    newy = np.zeros(new_shape)
    newx[nhalo:new_shape[0] -nhalo, nhalo:new_shape[1] -nhalo] = np.transpose(x)
    newy[nhalo:new_shape[0] -nhalo, nhalo:new_shape[1] -nhalo] = np.transpose(y)
    #print fname
    b = SimulationBlock(2, block_number=no)
    g1 = h5f.create_group(b.blockname)
    apply_group_attributes(g1, b)
    #block_dset_name = b.location_dataset("x0").base
    block_dset_name = "x0_B%d" %no
    dset = g1.create_dataset('%s' % (block_dset_name), data=newx)
    set_hdf5_metadata(dset, halos=[-nhalo, nhalo], npoints=[shape[0], shape[1]], block=b)
    
    #block_dset_name = b.location_dataset("x1").base
    block_dset_name = "x1_B%d" %no
    dset = g1.create_dataset('%s' % (block_dset_name), data=newy)
    set_hdf5_metadata(dset, halos=[-nhalo, nhalo], npoints=[shape[0], shape[1]], block=b)
    
    print "LX block :", abs(np.amin(x) - np.amax(x))
    print "Ly block:", abs(np.amin(y) - np.amax(y))
#h5f.close()
    
