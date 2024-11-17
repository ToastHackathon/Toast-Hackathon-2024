# -*- coding: utf-8 -*-
"""
Created on Sat Nov 16 14:38:30 2024

@author: itsfr
"""

import os

# Change to the directory where this script is located
os.chdir(os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import polyscope as ps
import polyscope.imgui as psim
import argparse, json
from PARKER_RIGID import *
from contact import *
from Collision import *

# Initialize polyscope, ground plane at zero
ps.init()
ps.set_automatically_compute_scene_extents(False)
ps.set_length_scale(10)
low = np.array((-2, -2., -2.))
high = np.array((2., 2., 2.))
ps.set_bounding_box(low, high)
ps.set_always_redraw(True)
ps.set_ground_plane_height_factor(-0.2,is_relative=True)

# process argumetns from command line and setup with selected json file
parser = argparse.ArgumentParser()
parser.add_argument("--file", type=str, default = "scenes/scene2.json")
args = parser.parse_args()
print(args.file)
data = json.load(open(args.file))
rigid_body_list = []
for body_desc in data['bodies']:
    rigid_body_list.append( RigidBody( body_desc ) )
sim_params = data['sim_parameters'][0]
gravity = np.array(sim_params.get('gravity',(0.0,0.0,0.0)))	# 3x1 vector
h = sim_params.get('dt',0.01)
mu = sim_params.get('mu',0)
substeps = sim_params.get('substeps',1)
is_running = sim_params.get('is_running',False)
check_collisions = sim_params.get('check_collisions',True)
show_contacts = sim_params.get('show_contacts',False)
stop_on_contact = sim_params.get('stop_on_contact',False)

# other setup before main loop
collision = Collision()
elapsed = 0
np.set_printoptions(formatter={'float_kind':"{:.2f}".format})

def main_display_loop():
    global is_running, check_collisions, show_contacts, stop_on_contact, substeps, gravity, mu, elapsed, h

    # MIKEY BAKER, MCGILL ID 261106790
    psim.TextUnformatted("Mikey Baker, McGill ID: 261106790")
    psim.TextUnformatted(sim_params['name'])
    if(psim.Button("Reset")):
        collision.reset()
        collision.update_display(show_contacts)
        for rb in rigid_body_list:
            rb.reset()
            rb.update_display()
            elapsed = 0
    do_step = False
    if(psim.Button("Step")):
        do_step = True
    _, is_running	= psim.Checkbox("Run", is_running)
    psim.TextUnformatted("Elapsed = " + str(elapsed))
    _, h = psim.SliderFloat("step size", h, v_min=0.001, v_max=0.1)
    _, substeps = psim.InputInt("substeps", substeps, step=1, step_fast=1)
    if substeps < 1:
        substeps = 1;
    _, check_collisions = psim.Checkbox("Check collisions", check_collisions )
    _, stop_on_contact = psim.Checkbox("Stop on contact", stop_on_contact)
    changed, show_contacts = psim.Checkbox("Show contacts", show_contacts)
    if changed:
        collision.update_display(show_contacts)
    _, gravity[1] = psim.SliderFloat("gravity y", gravity[1], v_min=-10, v_max=10)
    _, mu = psim.SliderFloat("friction", mu, v_min=0, v_max=5)
    showFirstBodyStats=False
    if (showFirstBodyStats):
        rb = rigid_body_list[0]
        psim.TextUnformatted(rb.name + ':')
        L = rb.mass*rb.omega
        kinetic = 0.5*rb.mass*np.dot(rb.v,rb.v) + 0.5*np.dot(rb.omega,L)
        p = rb.mass*rb.v
        potential = np.dot(rb.x[1],gravity.T)*rb.mass
        psim.TextUnformatted("Angular momentum = " + str(L) + "total = " + str(np.linalg.norm(L)))
        psim.TextUnformatted("Kinetic energy = " + str(kinetic))
        psim.TextUnformatted("Linear momentum = " + str(p)+ "total = " + str(np.linalg.norm(p)))
        psim.TextUnformatted("Grav. Potential Energy = " + str(potential))

    # show the energy and momentum of just the first body, only a test option
    if(psim.TreeNode("Energy and Momentum")):
        for rb in rigid_body_list:
            psim.TextUnformatted(rb.name)
            #TODO: compute and display the kinetic energy, potential energy, and linear and angular momentum of each body
            #psim.TextUnformatted(rb.name + ':')
            L = np.matmul(rb.J0,rb.omega)
            kinetic = 0.5*rb.mass*np.dot(rb.v,rb.v) + 0.5*np.dot(rb.omega,L)
            p = rb.mass*rb.v
            potential = np.dot(rb.x[1],gravity.T)*rb.mass
            psim.TextUnformatted("Angular momentum = " + str(L) + "total = " + str(np.linalg.norm(L)))
            psim.TextUnformatted("Kinetic energy = " + str(kinetic))
            psim.TextUnformatted("Linear momentum = " + str(p)+ "total = " + str(np.linalg.norm(p)))
            psim.TextUnformatted("Grav. Potential Energy = " + str(potential))




        psim.TreePop()

    if is_running or do_step:
        # ignore substeps for now
        for i in range(substeps):
            stepsize = h/substeps
            for rb in rigid_body_list:
                rb.zero_force_and_torque()
                rb.add_force( gravity * rb.mass )
                rb.step_vel( stepsize )
            if check_collisions and collision.check(rigid_body_list) :
                collision.process(rigid_body_list,mu,stepsize)
                if stop_on_contact:
                    is_running = False
            for rb in rigid_body_list:
                rb.step_pos( stepsize )
            elapsed += stepsize
        for rb in rigid_body_list:
            rb.update_display()
        if show_contacts:
            collision.update_display(show_contacts)

#Parameters that are seen by the user when opening polyscope
is_true1 = False
is_running = False
is_true2 = True
ui_bread_type = {"Whole wheat bread":325,
                   "White bread": 200,
                   "Brick bread": 2000,
                   "Banana bread": 500}  # in kg/m^3
ui_bread_type_selected = list(ui_bread_type.keys())[1]
ui_bread_density = ui_bread_type["Whole wheat bread"]
ui_spread_type = ["Butter", "Jam", "Peanut butter"]
ui_spread_type_selected = ui_spread_type[0]
ui_shape_type = ['Toast','Triangle toast', 'Bagel', 'Croissant']
ui_shape_type_selected = ui_shape_type[0]
ui_initial_position = (0.0, 0.0, 1.0)
ui_initial_velocity = (0.0,0.0,0.0)
ui_initial_angle= 0.0
ui_initial_angular_velocity= (0.0,0.0,0.0)
#ui_toastiness = 0
bread_colors = {'rare': np.array([219.0, 177.0, 118.0]) / 255,
                'med-rare': np.array([198.0, 134.0, 78.0]) / 255,
                'med': np.array([189.0, 97.0, 35.0]) / 255,
                'well-done': np.array([97.0, 27.0, 8.0]) / 255,
                'burnt': np.array([16.0, 18.0, 17.0]) / 255}


First = True
# Define our callback function, which Polyscope will repeatedly execute while running the UI.
# We can write any code we want here, but in particular it is an opportunity to create ImGui 
# interface elements and define a custom UI.
ps.set_build_gui(False)
def callback():

    # If we want to use local variables & assign to them in the UI code below, 
    # we need to mark them as nonlocal. This is because of how Python scoping 
    # rules work, not anything particular about Polyscope or ImGui.
    # Of course, you can also use any other kind of python variable as a controllable 
    # value in the UI, such as a value from a dictionary, or a class member. Just be 
    # sure to assign the result of the ImGui call to the value, as in the examples below.
    # 
    # If these variables are defined at the top level of a Python script file (i.e., not
    # inside any method), you will need to use the `global` keyword instead of `nonlocal`.
    global ui_bread_density, First,bread_colors,is_running,is_true1, is_true2, ui_initial_position, ui_bread_type_selected, ui_spread_type_selected, ui_shape_type_selected, ui_initial_velocity, ui_initial_angle, ui_initial_angular_velocity, ui_toastiness


    # == Settings

    psim.PushItemWidth(150)
    

    #  Show text in the UI

    psim.TextUnformatted("Pick your bread type, spread, and initial conditions")
    
    psim.Separator()

    

#bread type
    #bread type
    psim.PushItemWidth(200)
    changed = psim.BeginCombo("Pick one type of bread", ui_bread_type_selected)
    if changed:
        for type in ui_bread_type.keys():
            _, selected = psim.Selectable(type, ui_bread_type_selected==type)
            if selected:
                ui_bread_type_selected = type
                ui_bread_density = ui_bread_type[type]
        psim.EndCombo()
    psim.PopItemWidth()

#spread type 
    psim.PushItemWidth(200)
    changed = psim.BeginCombo("Pick one type of spread", ui_spread_type_selected)
    if changed:
        for val in ui_spread_type:
            _, selected = psim.Selectable(val, ui_spread_type_selected==val)
            if selected:
                ui_spread_type_selected = val
        psim.EndCombo()
    psim.PopItemWidth()

#shape of bread
    psim.PushItemWidth(200)
    changed = psim.BeginCombo("Pick one type of shape", ui_shape_type_selected)
    if changed:
        for val in ui_shape_type:
            _, selected = psim.Selectable(val, ui_shape_type_selected==val)
            if selected:
                ui_shape_type_selected = val
        psim.EndCombo()
    psim.PopItemWidth()
    #toastiness
    rb = rigid_body_list[0]

    if(psim.Button("Rare")):
        rb.ps_mesh.set_color(tuple(bread_colors['rare']))
    psim.SameLine()

    if(psim.Button("Medium Rare")):
        rb.ps_mesh.set_color(tuple(bread_colors['med-rare']))
    psim.SameLine()

    if(psim.Button("Medium")):
        rb.ps_mesh.set_color(tuple(bread_colors['med']))
    psim.SameLine()

    if(psim.Button("Well Toasted")):
        rb.ps_mesh.set_color(tuple(bread_colors['well-done']))
    psim.SameLine()

    if(psim.Button("Burnt")):
        rb.ps_mesh.set_color(tuple(bread_colors['burnt']))

    psim.TextUnformatted("Now, pick your initial position and velocity as well as the inclination angle")
    
    psim.Separator()
    



#initial position
    psim.PushItemWidth(300)
    changed, ui_initial_position = psim.InputFloat3('Initial Position', ui_initial_position)
#initial velocity 
    psim.PushItemWidth(300)
    changed, ui_initial_velocity = psim.InputFloat3('Initial Velocity', ui_initial_velocity)
#initial inclination angle 
    changed, ui_initial_angle = psim.SliderAngle("Inclination Angle", ui_initial_angle, 
                v_degrees_min=-90, v_degrees_max=90)
#initial angular velocity
    psim.PushItemWidth(300)
    changed, ui_initial_angular_velocity = psim.InputFloat3('Initial Angular Velocity', ui_initial_angular_velocity)

    if(First):
        for rb in rigid_body_list:
            rb.v = ui_initial_velocity
            rb.x = ui_initial_position
            rb.omega = ui_initial_angular_velocity
            th = ui_initial_angle
            rb.R = np.array([[1,0,0],[0,np.cos(th),-np.sin(th)],[0,np.sin(th),np.cos(th)]])
            rb.update_display()
            rb.mass = ui_bread_density*rb.Volume
# toastiness


#Run button
    changed, is_running=psim.Checkbox("Drop the Toast",is_running)
    if(is_running):
        # ignore substeps for now
        for rb in rigid_body_list:
            rb.elapsed += h

            #rb.zero_force_and_torque()
                #rb.add_force( gravity * rb.mass )
            rb.updateVelocity()
        #if check_collisions and collision.check(rigid_body_list) :
                #collision.process(rigid_body_list,mu,stepsize)
          #  i=i
                #if stop_on_contact:
                #    is_running = False
        for rb in rigid_body_list:
            rb.updatePosition()
            #elapsed += stepsize
        for rb in rigid_body_list:
            if rb.rad_crit>rb.x[1]:
                #print('calling resolve contact')
                bool = rb.resolveContact()
            if rb.run:
                rb.update_display()
        #if show_contacts:
        #    collision.update_display(show_contacts)
        First = False

ps.set_user_callback(callback)
ps.show()


