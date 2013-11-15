beh-zee-yay
===========
Names: Jason Ye, Amy Tsai

Platform: OSX 

Location: Amy Tsai is submitting onto the OSX machines.

Building the Project
--------------------
type 
'''
make
'''
into the command line

Command line arguments
----------------------
After making the project, o run the program type:
'''
./bezier <.bez file> <subdivision parameter> [-a]
'''

* The subdivision parameter specifies the __step size__ for uniform subdivision and the allowable __error__ for adaptive subdivision
* By default, the program runs with uniform subdivision, -a runs adaptive

Navigation
----------
Once the program is launched, the default view shows the entire scene the following commands can be used for navigation:

* Arrow keys __rotate__ the scene 
* SHIFT + arrow keys translate the scene
* '+' and '-' (keys without shift) zoom the scene in and out

Shading Options
---------------
The default view is smooth Gouraud shading. The following commands can be used to change the shading:

* 's' toggles between __smooth__ and __flat__ shading
* 'w' toggles between __wireframe__ and __filled__ polygons

Extra Features
--------------

Our program allows for the interactive chaingg of the input parameter by holding down ALT and pressing any of the arrow keys
* Up: Parameter multiplied by 2
* Down: Parameter multiplied by .5
* Right: Parameter multiplied by 1.5
* Left: Parameter multiplied by .66