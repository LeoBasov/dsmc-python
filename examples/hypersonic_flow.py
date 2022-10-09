import dsmc
import dsmc.writer as wrt
import dsmc.diagnostics as dia
import dsmc.mesh.mesh2d as msh2d
import dsmc.octree as oc
import time
import numpy as np

if __name__ == '__main__':
    # general parameters
    solver = dsmc.DSMC()
    domain = [(-3.0, 3.0), (-1.5, 1.5), (-0.025, 0.025)]
    obj = [(-0.25, 0.25), (-0.25, 0.25), (-0.5, 0.5)]
    dt = 1e-6
    w = 0.2 * 2.3e+15
    mass = 6.6422e-26
    T =  273.0
    n = 2.6e+19
    u = np.array([0.0, -3043.0, 0.0])
    niter = 1000
    niter_sample = 100
    
    # set up mesh2
    mesh = msh2d.Mesh2d()
    
    mesh.n_cells1 = 100
    mesh.n_cells2 = 50
    mesh.min1 = domain[0][0]
    mesh.min2 = domain[1][0]
    mesh.cell_size1 = 0.06
    mesh.cell_size2 = 0.06
    
    h = domain[2][1] - domain[2][0]
    
    # setup solver
    solver.set_domain(domain)
    solver.set_weight(w)
    solver.set_mass(mass)
    
    # set up boundary
    solver.set_bc_type("xmin", "inflow")
    solver.set_bc_type("xmax", "inflow")
    solver.set_bc_type("ymax", "inflow")
    
    solver.set_bc_type("ymin", "open")
    
    solver.set_bc_values("xmin", T, n, u)
    solver.set_bc_values("xmax", T, n, u)
    solver.set_bc_values("ymax", T, n, u)
    
    # set object
    solver.add_object(obj)
    
    # create particles
    solver.create_particles(domain, T, n, u)
    
    # trees
    tree_inner = oc.Octree()
    tree_outer = oc.Octree()
 
    inner_pos = np.array([[-0.25, -0.25, -0.025], [0.25, 0.25, 0.025]])
    outer_pos = np.array([[-3.0, -1.5, -0.025], [3.0, 1.5, 0.025]])
    
    tree_inner.build(inner_pos)
    tree_outer.build(outer_pos)
  
    wrt.write_buttom_leafs(tree_inner, "innter.vtu")
    wrt.write_buttom_leafs(tree_outer, "outer.vtu")
    
    # start timing
    start_time = time.time()
    
    for it in range(niter):
        print("iteration {:4}/{}".format(it + 1, niter), end="\r", flush=True)
        solver.advance(dt)
            
    for it in range(niter_sample):
        print("iteration {:4}/{}".format(it + 1, niter_sample), end="\r", flush=True)
        solver.advance(dt)
        
        if (it + 1)%10 == 0:
            mesh.sort(solver.particles.Pos)
            boxes, values = dia.analyse_2d(solver.particles.Pos, solver.particles.Vel, mesh, h)
            wrt.write_planar(boxes, values, "hypersonic_{}.vtu".format(it))
        
    wrt.write_buttom_leafs(solver.octree)
        
    print("")
    print("--- %s seconds ---" % (time.time() - start_time))
    print('done')