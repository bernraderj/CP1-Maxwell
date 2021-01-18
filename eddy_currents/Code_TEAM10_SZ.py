from dolfin import *
import numpy as np
from petsc4py import PETSc
from ufl import sign


set_log_level(15)

Omega = {'air':3, 'coil':1, 'conductor':2, 'boundary':4}

scale = 1e-3
mesh = Mesh("Mesh_TEAM10.xml")
domains = MeshFunction("size_t", mesh, 3, mesh.domains())
dx = Measure("dx", subdomain_data=domains)

File("domains.pvd") << domains

class CurrentExpression(UserExpression):
  def __init__(self, size, j = 1.0968, center = None, **kwargs): #scaling factor
    super().__init__()
    self._j = j
    self._center = center
    self._size = size

  def eval(self, value, x):
    # handle offset
    xc = x.copy()
    if self._center is not None:
      xc -= np.array(self._center)

    value[2] = 0.0
    # parametrize straigth lines
    size = self._size
    if -size[0]/2. <= xc[0] <= size[0]/2. and xc[1] > 0.:
        value[0] = -1.0
        value[1] =  0.0
    if -size[0]/2. <= xc[0] <= size[0]/2. and xc[1] < 0.:
        value[0] =  1.0
        value[1] =  0.0
    if xc[0] > 0. and -size[1]/2. <= xc[1] <= size[1]/2.:
        value[0] =  0.0
        value[1] =  1.0
    if xc[0] < 0. and -size[1]/2. <= xc[1] <= size[1]/2.:
        value[0] =  0.0
        value[1] = -1.0

    # parametrize curves
    if xc[0] > size[0]/2. and xc[1] > size[1]/2.:
        rx, ry = xc[0]-size[0]/2., xc[1]-size[1]/2.
        r = np.sqrt(rx**2 + ry**2)
        value[0] = -ry/r
        value[1] =  rx/r
    if xc[0] < -size[0]/2. and xc[1] > size[1]/2.:
        rx, ry = xc[0]+size[0]/2., xc[1]-size[1]/2.
        r = np.sqrt(rx**2 + ry**2)
        value[0] = -ry/r
        value[1] =  rx/r
    if xc[0] > size[0]/2. and xc[1] < -size[1]/2.:
        rx, ry = xc[0]-size[0]/2., xc[1]+size[1]/2.
        r = np.sqrt(rx**2 + ry**2)
        value[0] = -ry/r
        value[1] =  rx/r
    if xc[0] < -size[0]/2. and xc[1] < -size[1]/2.:
        rx, ry = xc[0]+size[0]/2., xc[1]+size[1]/2.
        r = np.sqrt(rx**2 + ry**2)
        value[0] = -ry/r
        value[1] =  rx/r

    value *= self._j

  def value_shape(self):
    return (3,)

V = FunctionSpace(mesh, "N1curl", 1)
VV = VectorFunctionSpace(mesh, "CG", 1)
u = TrialFunction(V)
v = TestFunction(V)

#center=np.array([194,100,0]),\

j_expr = CurrentExpression(size=np.array([100,100,0]), center=np.array([0, 0, 0]),\
    j=5.64*133/(0.1*0.025), degree=1) #2742./(0.1*0.025)
j_expr = interpolate(j_expr, VectorFunctionSpace(mesh, "CG", 1))

a = inner(curl(u),curl(v))*dx
L = Constant(scale)*inner(j_expr,curl(v))*dx(Omega['coil'])

bc = DirichletBC(V, Constant((0.,0.,0.)), DomainBoundary())

T0 = Function(V)
A, rhs = assemble_system(a, L, bc)

# Create PETSc Krylov solver (from petsc4py)
ksp = PETSc.KSP()
ksp.create(PETSc.COMM_WORLD)

# Set the Krylov solver type and set tolerances
ksp.setType("cg")
ksp.setTolerances(rtol=1.0e-14, atol=1.0e-14, divtol=1.0e10, max_it=300)

# Get the preconditioner and set type (HYPRE AMS)
pc = ksp.getPC()
pc.setType("hypre")
pc.setHYPREType("ams")

# Build discrete gradient
P1 = FunctionSpace(mesh, "CG", 1)
G = DiscreteOperators.build_gradient(V, P1)

# Attach discrete gradient to preconditioner
pc.setHYPREDiscreteGradient(as_backend_type(G).mat())

# Build constants basis for the Nedelec space
constants = [Function(V) for i in range(3)]
for i, c in enumerate(constants):
    direction = [1.0 if i == j else 0.0 for j in range(3)]
    c.interpolate(Constant(direction))

# Inform preconditioner of constants in the Nedelec space
cvecs = [as_backend_type(constant.vector()).vec() for constant in constants]
pc.setHYPRESetEdgeConstantVectors(cvecs[0], cvecs[1], cvecs[2])

# We are dealing with a zero conductivity problem (no mass term), so
# we need to tell the preconditioner
pc.setHYPRESetBetaPoissonMatrix(None)

# Set operator for the linear solver
ksp.setOperators(as_backend_type(A).mat())

# Set options prefix
ksp.setOptionsPrefix("eddy_")

# Turn on monitoring of residual
opts = PETSc.Options()
opts.setValue("-eddy_ksp_monitor_true_residual", None)

# Solve eddy currents equation (using potential T0)
ksp.setFromOptions()
ksp.solve(as_backend_type(rhs).vec(), as_backend_type(T0.vector()).vec())

# Show linear solver details
#ksp.view()

print('"cg", "ams", norm(A*T0.vector()-rhs) / norm(rhs), norm(A*T0.vector()-rhs) / norm(rhs)')
File("data/T0.pvd") << project(T0, VV, solver_type='cg', preconditioner_type='amg' )

#Now the eddy current problem should be solved with the calculated T0
f = 1e-1
omega = 2.*np.pi*f

Hcurl = FiniteElement("N1curl", tetrahedron, 1)
CG = FiniteElement("CG", tetrahedron, 1)
V = FunctionSpace(mesh, MixedElement((Hcurl, Hcurl)))
VV0 = VectorFunctionSpace(mesh, "CG", 1)
V0 = FunctionSpace(mesh, "DG", 0)

u_r, u_i = TrialFunctions(V)
v_r, v_i = TestFunctions(V)

mu0 = 4.*np.pi*1e-7
mus = {1: mu0, 2: 2e3*mu0, 3: 1*mu0}
mu = Function(V0)
mu.vector()[:] = np.array([mus[i] for i in domains.array()])

sigmas = {1: 1e1, 2: 7.505e6, 3: 1e1}
sigma = Function(V0)
sigma.vector()[:] = np.array([sigmas[i] for i in domains.array()])

a =  1/mu*inner(curl(u_r), curl(v_r))*dx
a -= 1/mu*inner(curl(u_i), curl(v_i))*dx

a -= Constant(scale**2*omega)*sigma*inner(u_i, v_r)*dx
a -= Constant(scale**2*omega)*sigma*inner(u_r, v_i)*dx

L = inner(T0, curl(v_r))*dx

bcs = [ DirichletBC(V.sub(0), Constant((0., 0., 0.)), DomainBoundary()),
        DirichletBC(V.sub(1), Constant((0., 0., 0.)), DomainBoundary())]

A = Function(V)
M, rhs = assemble_system(a, L, bcs)

p = 1/mu*inner(curl(u_r), curl(v_r))*dx
p += 1/mu*inner(curl(u_i), curl(v_i))*dx
p += Constant(scale**2)*Constant(omega)*sigma*inner(u_r, v_r)*dx
p += Constant(scale**2)*Constant(omega)*sigma*inner(u_i, v_i)*dx
P, _ = assemble_system(p, L, bcs)

ksp = PETSc.KSP()
ksp.create(PETSc.COMM_WORLD)

ksp.setType("fgmres")
ksp.setTolerances(rtol = 1.0e-15, atol = 1.0e-15, max_it = 20)

P1 = FunctionSpace(mesh, "CG", 1)
G = as_backend_type(DiscreteOperators.build_gradient(FunctionSpace(mesh, "N1curl", 1), P1)).mat()

V1 = FunctionSpace(mesh, "N1curl", 1)
constants = [Function(V1) for i in range(3)]
for i, c in enumerate(constants):
    direction = [1.0 if i == j else 0.0 for j in range(3)]
    c.interpolate(Constant(direction))

pc = ksp.getPC()
pc.setType("fieldsplit")
pc.setFieldSplitType(0)
is0 = PETSc.IS().createGeneral(V.sub(0).dofmap().dofs())
is1 = PETSc.IS().createGeneral(V.sub(1).dofmap().dofs())

pc.setFieldSplitIS(('A_r', is0), ('A_i', is1))
for subksp in pc.getFieldSplitSubKSP():
    subksp.setType("cg")
    subksp.setTolerances(rtol = 1.0e-14, atol=1.0e-14, max_it = 6)
    subpc = subksp.getPC()
    subpc.setType("hypre")
    subpc.setHYPREType("ams")
    subpc.setHYPREDiscreteGradient(G)
    cvecs = [as_backend_type(constant.vector()).vec() for constant in constants]
    subpc.setHYPRESetEdgeConstantVectors(cvecs[0], cvecs[1], cvecs[2])

ksp.setOperators(as_backend_type(M).mat(), as_backend_type(P).mat())

ksp.setOptionsPrefix("eddy_")

opts = PETSc.Options()
opts.setValue("-eddy_ksp_monitor_true_residual", None)

ksp.setFromOptions()
ksp.solve(as_backend_type(rhs).vec(), as_backend_type(A.vector()).vec())

Ar_plate, Ai_plate = A.split()
Br = project(curl(Ar_plate), VV0, solver_type = 'cg', preconditioner_type = 'amg')
Bi = project(curl(Ai_plate), VV0, solver_type = 'cg', preconditioner_type = 'amg')
Hr = project((1/mu)*curl(Ar_plate), VV0, solver_type = 'cg', preconditioner_type = 'amg')
Hi = project((1/mu)*curl(Ai_plate), VV0, solver_type = 'cg', preconditioner_type = 'amg')

#mesh_plate = SubMesh(mesh, Omega['conductor'])
#Ai_plate = interpolate(A_i, FunctionSpace(mesh_plate, "N1curl", 1))
#Ar_plate = interpolate(A_r, FunctionSpace(mesh_plate, "N1curl", 1))

Jr = project(Constant(scale*omega)*sigma*Ai_plate,\
    VectorFunctionSpace(mesh, "CG", 1), solver_type='cg',\
    preconditioner_type='ilu')

#Jr1 = project(Constant((0, 0, 1e3)) + Jr, VectorFunctionSpace(mesh_plate, "CG", 1), solver_type='cg',\
#    preconditioner_type='ilu')

Ji = project(-Constant(scale*omega)*sigma*Ar_plate,\
    VectorFunctionSpace(mesh, "CG", 1), solver_type='cg',\
    preconditioner_type='ilu')

J = project(sigma*Constant(omega*scale)*sign(Ai_plate[1])*sqrt((Ai_plate[1])**2 + (Ar_plate[1])**2),\
    FunctionSpace(mesh, "CG", 1),\
    solver_type='cg', preconditioner_type = 'ilu')

File("data/AAstar_Ar.pvd") << Ar_plate
File("data/AAstar_Ai.pvd") << Ai_plate

File("data/AAstar_Jr.pvd") << Jr
#File("data/AAstar_Jr1_50Hz.pvd") << Jr1
File("data/AAstar_ji.pvd") << Ji
File("data/AAstar_J.pvd") << J
File("data/AAstar_Br.pvd") << Br
File("data/AAstar_Bi.pvd") << Bi
File("data/AAstar_Hr.pvd") << Hr
File("data/AAstar_Hi.pvd") << Hi
