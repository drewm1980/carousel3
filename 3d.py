#!/usr/bin/env python
import numpy
from numpy import linspace, pi, array, cos, sin, zeros, array, sqrt, sign
from numpy.linalg import solve, norm
#import casadi
#from casadi import ssym, vertcat

# This python script solves a 3d version of the carousel3
# tether cost optimization problem.
# For the 3D case, we need two rings of anchors.
# Otherwise the hitch will move off the vertical axis, screwing up the geometry.

# Much of the geometry can be parameterized by circles at ground level.
# These include, in order from large to small,
# * Outer Anchor Circle   -> r_outer_anchor
# * Upper guy Force Cone, -> r_force_cone
# * Inner Anchor Circle   -> r_inner_anchor

# The most import design criterion
r_boom = 10. # Radius of boom end circle of travel in m
h_boom = 8. # height of boom end circle of travel in m

anchor_count = 4 # Number of anchors per anchor ring
r_force_cone = 10. # in m, Indirectly influences hitch height
eps=.001
ratio_outer_anchor = 1.0-eps # ratio of r_outer_anchor_incircle/r_force_cone.  0 to 1-eps
ratio_inner_anchor = 1.0-eps # ratio of r_force_cone/r_inner_anchor.  0 to 1-eps
f_tether = 1. # N, radial direction
tether_cost_per_volume = 1. # euros per unit of volume
tether_breaking_strength = 1. # units of pressure

def normalize(v):
    return v/norm(v)

# Solution for circle-line intersection lifted from Mathworld
def intersect_circle_line(x1,y1,x2,y2,r):
    dx = x2-x1
    dy = y2-y1
    dr = sqrt(dx*dx + dy*dy)
    D = x1*y2 - x2*y1
    delta = r*r*dr*dr - D*D
    if delta <=0:
        raise Exception('Description has wrong sign; no intersection')
    xa = (D*dy + sign(dy)*dx*sqrt(delta))/(dr*dr)
    xb = (D*dy - sign(dy)*dx*sqrt(delta))/(dr*dr)
    ya = (-D*dx + abs(dy)*sqrt(delta))/dr*dr
    yb = (-D*dx - abs(dy)*sqrt(delta))/dr*dr
    return xa,xb,ya,yb

def analyze_carousel(r_boom, 
                     h_boom, 
                     anchor_count, 
                     r_force_cone,
                     ratio_outer_anchor,
                     ratio_inner_anchor,
                     f_tether,
                     tether_cost_per_volume,
                     tether_breaking_strength):
    r_anchor_outer = r_force_cone / (sin(pi/anchor_count)*ratio_outer_anchor)
    r_anchor_inner = r_force_cone * ratio_inner_anchor
    l_boom = norm([r_boom,h_boom]) # length of the boom in m

    h_hitch = h_boom/(r_boom+r_force_cone)*r_force_cone

    # convention: guy directions always point out of the hitch

    # Force ballance at the end of the boom, in vertical the 2D plane
    # passing through the boom and the active upper guy line
    guy0_dir = normalize(array([r_boom, h_boom-h_hitch])) # in boom plane
    boom_dir = normalize(array([r_boom,h_boom]))
    tether_dir = normalize(array([1.0,0.0]))
    # f_tether*tether_dir + f_guy0*(-guy0_dir) + f_boom*boom_dir = 0
    dirs = zeros([2,2])
    dirs[:,0] = -guy0_dir
    dirs[:,1] = boom_dir
    guy0_and_boom_force = solve(dirs, -tether_dir*f_tether)
    f_guy0 = guy0_and_boom_force[0] # This a tension
    f_boom = guy0_and_boom_force[1] # This is a compression
    assert(f_guy0>-1e-9)
    assert(f_boom>-1e-9)

    # Horizontal angle between consecutive outer and inner anchors, in radians
    a_anchor = pi/anchor_count

    # Some geometric constructions to get the transition angle between two
    # cells in closed form
    slope = r_anchor_inner*sin(a_anchor)/(r_anchor_outer - r_anchor_inner*cos(a_anchor))
    # r_force_cone * sin(theta_trans) = (r_anchor - r_force_cone*cos(theta_trans)) * slope
    # sin(theta_trans) = (r_anchor_outer - r_force_cone * cos(theta_trans))*slope/r_force_cone
    # sin(theta_trans) = r_anchor_outer*slope/r_force_cone - slope*cos(theta_trans)
    # sin(theta_trans)+slope*cos(theta_trans) = r_anchor_outer*slope/r_force_cone

    # Working here!
    # Angle of rotation is counter-clockwise from above.
    # at theta = 0, boom starts aligned with negative x axis
    # There is an outer anchor aligned with positive x axis
    xa,xb,ya,yb = intersect_circle_line(r_anchor_inner*cos(a_anchor),
                          r_anchor_inner*sin(a_anchor),
                          r_anchor_outer,
                          0,
                          r_force_cone)
    print xa,xb,ya,yb
    
    #max_guy_forces = zeros((2,1))
    #for delta in linspace(0.0,anchorAngle):
        ## Force ballance at centerpoint
        #guy0_dir = array([cos(delta),sin(delta)])
        #a1 = pi
        #a2 = pi+anchorAngle
        #guy1_dir = array([cos(a1),sin(a1)])
        #guy2_dir = array([cos(a2),sin(a2)])
        #guy12_dirs = zeros([2,2])
        #guy12_dirs[:,0] = guy1_dir
        #guy12_dirs[:,1] = guy2_dir
        ## guy1_dir*f_guy1 + guy2_dir*f_guy2 + guy0_dir*f_guy0 = 0
        #guy_forces = solve(guy12_dirs,-guy0_dir*f_guy0)
        #assert(guy_forces[0]>-1e-9)
        #assert(guy_forces[1]>-1e-9)
        #for i in range(2):
            #if guy_forces[i] > max_guy_forces[i]:
                #max_guy_forces[i] = guy_forces[i]
                ## Compute system cost
                ##boom_cost = l_boom*f_boom
                #boom_cost = 0 # Just to take it out of the budget for now.
                #l_guy0 = r_boom # Approximation!
                #guy0_cost = l_guy0*f_guy0*tether_cost_per_volume*tether_breaking_strength
                #l_guy12 = r_ring # Approximation!
                #guy12_cost = l_guy12 * max_guy_forces[0] * tether_cost_per_volume * tether_breaking_strength
                ## Note, this assumes we're "splitting" the web twice to allow the boom
                ## to pass throught the web.
                #if anchor_count == 3:
                    #num_splits = 2
                #else:
                    #num_splits = 2
                    #system_cost = boom_cost + num_splits*guy0_cost + (anchor_count+num_splits)*guy12_cost
        #return anchor_count, r_ring, num_splits, float(system_cost)


if __name__=='__main__':
    for anchor_count in range(4,10,2):
        for r_outer_anchor in linspace(0.5*r_boom + pi*1e-6,1.2*r_boom):
            for r_inner_anchor in linspace(0,r_outer_anchor):
                analyze_carousel(r_boom, 
                                 h_boom, 
                                 anchor_count, 
                                 r_force_cone,
                                 ratio_outer_anchor,
                                 ratio_inner_anchor,
                                 f_tether,
                                 tether_cost_per_volume,
                                 tether_breaking_strength)

