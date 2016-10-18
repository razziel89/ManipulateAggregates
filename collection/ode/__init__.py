"""
function definitions for advacing out little world using ode
This is just a collection of usefull wrapper functions for ode
"""
#This file is part of ManipulateAggregates.
#
#Copyright (C) 2016 by Torsten Sachse
#
#ManipulateAggregates is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#ManipulateAggregates is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with ManipulateAggregates.  If not, see <http://www.gnu.org/licenses/>.
import ode
def transfer_positions(ode_bodies):
    """
    This functions returns the positions of the ode bodies "ode_bodies".
    This allows for easily updating what's to be displayed.
    """
    return [body.getPosition() for body in ode_bodies]

def ode_create_sphere(world,space,density,p,r):
    """
    Create an ode sphere at the coordinates "p" (3-dim vector) in the space "space"
    in the world "world" with the density "density" and the radius "r".
    """
    #create a body that obeys the worlds rules
    body = ode.Body(world)
    #initialize empty mass object
    M = ode.Mass()
    #set mass equal to that of a sphere with the given density and radius
    M.setSphere(density,r)
    #assign mass to object
    body.setMass(M)
    #collision detection is handled in the space using different collision geometries
    #this has nothing to do with the world the objects reside in
    #an object can have multiple collision geometries assigned to it but exist only once
    #in the world
    #create geometry that obeys the rules of the space
    geom = ode.GeomSphere(space=space,radius=r)
    #assign the geometry to the body
    geom.setBody(body)
    #position the body where it's supposed to be
    body.setPosition(p)
    return body, geom

def ode_near_callback(args, geom1, geom2):
    """
    This function checks if the given geoms do collide and
    creates contact joints if they do.
    """
    #assign variables hidden in "args"
    world,contactgroup = args
    # Check if the objects do collide
    contacts = ode.collide(geom1, geom2)
    # Create contact joints
    for c in contacts:
        #restitution parameter: basically how much speed is not lost upon impact (0: completely inelastic, 1: completely elastic)
        c.setBounce(0.1)
        #this is the Coulomb friction parameter. This is approx. steel on steel
        c.setMu(0.5)
        #create a contact joint for each contact and assign the corresponding bodies to it
        j = ode.ContactJoint(world, contactgroup, c)
        j.attach(geom1.getBody(), geom2.getBody())

def ode_next_step(world,space,contactgroup,dt,near_callback):
    """
    Advance the ode world by one timestep.
    world: the world where the simulation resides
    space: the space where the simulation resides
    contactgroup: the contactgroup used to find out whether objects collide
    dt: timestep for the simulation
    near_callback: function name used to find out whether things collide
    """
    # Detect collisions and create contact joints
    space.collide((world,contactgroup), near_callback)
    # Simulation step
    world.step(dt)
    # Remove all contact joints
    contactgroup.empty()
    return
    

