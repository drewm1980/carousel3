#!/usr/bin/env python
import numpy
from numpy import linspace, pi, array, cos, sin, zeros, array
from numpy.linalg import solve, norm
#import casadi
#from casadi import ssym, vertcat

# This python script solves a 2d version of the carousel3
# cost optimization problem.

# Simplifying Assumption
# When multiple guy lines are geometrically active,
# the lines that are most aligned with the tether will
# take up the full load.  In reality, the load may be 
# more distributed, so we are conservatively assuming
# worst-case loading.

r_boom = 10 # m
h_boom = 8 # m
l_boom = norm([r_boom,h_boom])
f_tether = 1 # N, radial direction

tether_cost_per_volume = 1 # euros per unit of volume
tether_breaking_strength = 1 # units of pressure

# For 2D case, single ring of anchors is sufficient.
for anchorCount in range(3,5):
    for r_ring in linspace(0.5*r_boom + pi*1e-6,1.2*r_boom):
        # Force ballance in vertical 2D plane guy attachment point
        boom_dir = array([r_boom,h_boom])
        boom_dir = boom_dir/norm(boom_dir)
        tether_dir = array([1.0,0.0])
        guy0_dir = array([-(r_boom+r_ring), -h_boom])
        guy0_dir = guy0_dir/norm(guy0_dir)
        # guy0_dir*f_guy0 + boom_dir*f_boom + tether_dir*f_tether = 0
        dirs = zeros([2,2])
        dirs[:,0] = guy0_dir
        dirs[:,1] = boom_dir
        guy0_and_boom_force = solve(dirs, -tether_dir*f_tether)
        f_guy0 = guy0_and_boom_force[0]
        f_boom = guy0_and_boom_force[1]
        assert(f_guy0>-1e-9)
        assert(f_boom>-1e-9)
        # From here out, do 2D analysis/approximation in horizontal plane
        anchorAngle = 2*pi/anchorCount # in horizontal plane
        max_guy_forces = zeros((2,1))
        for delta in linspace(0.0,anchorAngle):
            # Force ballance at centerpoint
            guy0_dir = array([cos(delta),sin(delta)])
            a1 = pi
            a2 = pi+anchorAngle
            guy1_dir = array([cos(a1),sin(a1)])
            guy2_dir = array([cos(a2),sin(a2)])
            guy12_dirs = zeros([2,2])
            guy12_dirs[:,0] = guy1_dir
            guy12_dirs[:,1] = guy2_dir
            # guy1_dir*f_guy1 + guy2_dir*f_guy2 + guy0_dir*f_guy0 = 0
            guy_forces = solve(guy12_dirs,-guy0_dir*f_guy0)
            assert(guy_forces[0]>-1e-9)
            assert(guy_forces[1]>-1e-9)
            for i in range(2):
                if guy_forces[i] > max_guy_forces[i]:
                    max_guy_forces[i] = guy_forces[i]
        # Compute system cost
        #boom_cost = l_boom*f_boom
        boom_cost = 0 # Just to take it out of the budget for now.
        l_guy0 = r_boom # Approximation!
        guy0_cost = l_guy0*f_guy0*tether_cost_per_volume*tether_breaking_strength
        l_guy12 = r_ring # Approximation!
        guy12_cost = l_guy12 * max_guy_forces[0] * tether_cost_per_volume * tether_breaking_strength
        # Note, this assumes we're "splitting" the web twice to allow the boom
        # to pass throught the web.
        if anchorCount == 3:
            num_splits = 2
        else:
            num_splits = 2
        system_cost = boom_cost + num_splits*guy0_cost + (anchorCount+num_splits)*guy12_cost
        print anchorCount, r_ring, num_splits, float(system_cost)


