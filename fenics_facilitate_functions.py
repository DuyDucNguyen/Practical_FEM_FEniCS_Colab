#from __future__ import print_function
#from fenics import *
#from facilitate_functions import myprint_green


def mesh_info(mesh):
	mesh.init()
	str = '\n \
	Number of vertices: {}\n \
	Number of edges:    {}\n \
	Number of faces:    {} \n \
	Number of facets:   {} \n \
	Number of cells:    {}'
	print( str.format(mesh.num_vertices(), mesh.num_edges(), mesh.num_faces(), mesh.num_facets(), mesh.num_cells()) )


def mpi_rank_size():
	#dolfin-version 2018.1.0
	comm = MPI.comm_world
	mpi_rank = MPI.rank(comm) # rank of current processor
	mpi_size = MPI.size(comm) # total processors called
	print("Processor {} is answering, of total {} processors".format(mpi_rank,mpi_size))
	return mpi_rank, mpi_size

mpi_rank, mpi_size = mpi_rank_size()


def msh2xml(input_dir, mesh_name):
	import os
	head_name, tail_name = mesh_name.split('.')
	myprint_green('>>> Converting mesh %s.msh into xml <<<' %(head_name))
	os.system("dolfin-convert %s.msh %s.xml" %(head_name, head_name))
	return


def xml2h5(filename = 'mesh.xml'):
	''' Description: Convert mesh written in XML format into H5 format
	    		 Write mesh tags (boundaries, domains) into H5 file
	    Input:  XML mesh
	    Output: H5 mesh
	'''
	head_name, tail_name = filename.split('.')
	mesh = Mesh(filename)
	cell_markers = MeshFunction("size_t", mesh, head_name + "_physical_region.xml")
	facet_markers = MeshFunction("size_t", mesh, head_name + "_facet_region.xml")
	hdf = HDF5File(mesh.mpi_comm(), head_name + ".h5", "w")
	hdf.write(mesh, "/mesh")
	hdf.write(cell_markers, "/cell_markers")
	hdf.write(facet_markers, "/facet_markers")
	#hdf.write(u, '/u')
	#print('xml2h5:', head_name + ".h5")
	return head_name + ".h5"





def read_h5(filename = 'file.h5'):
	'''Description: Convert mesh written in XML format into H5 format
	   Write mesh tags (boundaries, domains) into H5 file
	   Input:  XML mesh
	   Output: H5 mesh
	'''
	mesh = Mesh()
	hdf = HDF5File(mesh.mpi_comm(), filename, "r")
	hdf.read(mesh, "/mesh", False)
	dim = mesh.topology().dim()
	cell_markers = MeshFunction("size_t", mesh, dim)
	hdf.read(cell_markers, "/cell_markers")
	facet_markers = MeshFunction("size_t", mesh, dim-1)
	hdf.read(facet_markers, "/facet_markers")
	return mesh, cell_markers, facet_markers



def paraview_check_h5mesh(read_h5, wd = 'wd', mesh_name = 'mesh.h5', plot_option = True):
	'''Description: Read H5 mesh file
	   Input: H5 mesh file
	   Output: facet markers and cell markers in PVD format
	'''
	import os
	mesh, cell_markers, facet_markers = read_h5(wd+mesh_name)
	wd = wd + 'paraview_check_h5mesh/'
	if os.path.isdir(wd) == False:
		os.mkdir(wd)
	File(wd + "facet_markers.pvd") << facet_markers
	File(wd + "cell_markers.pvd") << cell_markers

	return mesh, cell_markers, facet_markers



def paraview_check_xml_mesh(
	input_dir = '',
	output_dir = '',
	mesh_filename = ''
	):
	"""
		DESCRIPTION:
		INPUT:
		OUTPUT:
	"""

	from fenics import Mesh, MeshFunction, File
	wd0 = output_dir + "Check_Geometry/"
	if os.path.isdir(wd0) == False:
		os.mkdir(wd0)

	mesh = Mesh(input_dir + mesh_filename + ".xml")

	# Use MeshFunction to read GMSH-tags of Boundaries
	facet_markers = MeshFunction("size_t", mesh, input_dir + mesh_filename + "_facet_region.xml")

	# Use MeshFunction to read GMSH-tags of sub_domains
	cell_markers  = MeshFunction("size_t", mesh, input_dir + mesh_filename +  "_physical_region.xml")

	# Save facet_markers
	File(wd0 + "facet_markers.pvd") << facet_markers
	File(wd0 + "cell_markers.pvd") << cell_markers

	return




def build_h5_mesh(**parameters):
	'''
		DESCRIPTION:
		INPUT:
		OUTPUT:
	'''
	if mpi_rank == 0:
		mesh_name = parameters['mesh_name']
		from create_mesh import create_gmsh_mesh
		head_name, tail_name = mesh_name.split('.')
		# Convert geo to msh
		#create_gmsh_mesh(**parameters)

		# Convert msh to xml
		msh2xml('', mesh_name)

		# Convert xml to h5
		h5mesh_name = xml2h5(filename = head_name + '.xml')
	return h5mesh_name


def extract_dofs_coor(mesh, degree):
	# Extract mesh coordinates
	mesh_coor = FunctionSpace(mesh, "Lagrange", degree).tabulate_dof_coordinates()
	mesh_coor = mesh_coor.reshape((-1, mesh.geometry().dim()))

	#mesh_coor_x = mesh_coor[:,0]
	#mesh_coor_y = mesh_coor[:,1]
	#mesh_coor_z = mesh_coor[:,2]
	return mesh_coor





def log_level_options():
	print("""
	CRITICAL  = 50, // errors that may lead to data corruption and suchlike
	ERROR     = 40, // things that go boom
	WARNING   = 30, // things that may go boom later
	INFO      = 20, // information of general interest
	PROGRESS  = 16, // what's happening (broadly)
	TRACE     = 13, // what's happening (in detail)
	""")

	#set_log_level(LogLevel.CRITICAL)
	#set_log_level(LogLevel.ERROR)
	#set_log_level(LogLevel.WARNING)
	#set_log_level(LogLevel.INFO)
	#set_log_level(LogLevel.PROGRESS)
	#set_log_level(LogLevel.TRACE)
	#set_log_level(True)

	return
