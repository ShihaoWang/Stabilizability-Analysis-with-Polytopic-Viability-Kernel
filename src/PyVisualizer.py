import sys, os, time
from klampt import *
from klampt import vis
from klampt.vis.glrobotprogram import GLSimulationPlugin
from klampt.model.trajectory import Trajectory
from scipy.interpolate import interp1d
import ipdb
import copy
from scipy.spatial import ConvexHull
import draw_hull
from OpenGL.GL import *
import math
import numpy as np

# This file contains functions related to

ExpName = "/home/motion/Desktop/Stabilizability-Analysis-with-Polytopic-Viability-Kernel/build/"
# ExpName = "/home/motion/Desktop/Scenarios/1/"
# ExpName = "/home/motion/Desktop/Scenarios/2/"
# ExpName = "/home/motion/Desktop/Scenarios/4/"

mode_no = 1;

class MyGLPlugin(vis.GLPluginInterface):
    def __init__(self, world):
        vis.GLPluginInterface.__init__(self)
        self.world = world
        self.quit = False
        self.starp = False
        self.mode_no = 1

    def mousefunc(self, button, state, x, y):
        print("mouse",button,state,x,y)
        if button==2:
            if state==0:
                print("Click list...",[o.getName() for o in self.click_world(x,y)])
            return True
        return False

    def motionfunc(self, x, y, dx, dy):
        return False

    def keyboardfunc(self, c, x, y):
        print("Pressed",c)
        if c == '1':
            self.mode_no = 1
            print "change mode no to 1"
            return True
        if c == '2':
            self.mode_no = 2
            print "change mode no to 2"
            return True
        if c == '3':
            self.mode_no = 3
            print "change mode no to 3"
            return True
        if c =='4':
            self.mode_no = 4
            print "change mode no to 4"
            return True
        return True

    def click_world(self, x, y):
        """Helper: returns a list of world objects sorted in order of
        increasing distance."""
        #get the viewport ray
        (s, d) = self.click_ray(x, y)

        #run the collision tests
        collided = []
        for g in self.collider.geomList:
            (hit, pt) = g[1].rayCast(s, d)
            if hit:
                dist = vectorops.dot(vectorops.sub(pt, s), d)
                collided.append((dist,g[0]))
        return [g[1] for g in sorted(collided)]

def my_draw_hull(h):
    glEnable(GL_LIGHTING)
    glEnable(GL_BLEND)
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA)
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,[1.0,0.25,0.5,0.5])
    draw_hull.draw_hull(h)

def String_List_to_Number_List(str_list, string_name):
    # This function is used to convert the string list to a certain type of list
    if string_name =="float":
        res_list = [float(i) for i in str_list]
    else:
        res_list = [int(i) for i in str_list]
    return res_list

def Configuration_Loader_fn(Config_Name):
    # This function is only used to load in the initial configuraiton
    # The initial file will be in the .config format
    with open(Config_Name,'r') as robot_angle_file:
        robotstate_angle_i = robot_angle_file.readlines()
    config_temp = [x.replace('\t',' ') for x in robotstate_angle_i]
    config_temp = [x.replace('\n','') for x in config_temp]
    config_temp = [float(i) for i in config_temp[0].split()]

    DOF = int(config_temp[0])
    # Config_Init = np.array(config_temp[1:])
    Config_Init = config_temp[1:]
    return DOF, Config_Init

def Read_Txt_fn(file_name):
    Empty_List = []
    with open(file_name) as Txt_File:
        Txt_File_Str = Txt_File.read().splitlines()
        for Txt_File_Str_i in Txt_File_Str:
            Txt_File_Str_i.replace("\'","")
            Empty_List.append(float(Txt_File_Str_i))
    return Empty_List

def State_Loader_fn(*args):
    if len(args) == 2:
        # In this case, the robot is only given the configuration file
        config_file_path = args[1] + args[0]
        DOF, Config_Init = Configuration_Loader_fn(config_file_path)
        # Then the Velocty_Init is set to be a zero value list
        Velocity_Init = []
        for i in range(0,DOF):
            Velocity_Init.append(0)
    else:
        if len(args) == 3:
            config_file_path = args[2] + args[0]
            velocity_file_path = args[2] + args[1]
            Config_Init = Read_Txt_fn(config_file_path)
            Velocity_Init = Read_Txt_fn(velocity_file_path)
            DOF = len(Config_Init)
        else:
            raise RuntimeError("Input name should be either one config file or two txt files!")
    return DOF, Config_Init, Velocity_Init

def Traj_Loader(robot_traj_path):
    with open(robot_traj_path,'r') as robot_traj_file:
        robot_traj_tot = robot_traj_file.readlines()
    return robot_traj_tot

def Traj_Loader_fn(*args):
    # This function is used to read in the traj file for visualization
    # However, this time we are going to load all the trajectories for visualization.
    # PVKRB, PVKCP, PVKHJB, ZSC, OE, CP,ZMP

    robot_traj_path = ""
    PVKRB_traj_path = ""
    PVKCP_traj_path = ""
    PVKHJB_traj_path = ""
    ZSC_traj_path = ""
    OE_traj_path = ""
    CP_traj_path = ""
    ZMP_traj_path = ""

    if len(args)==1:
        robot_traj_path = ExpName + "Data/stateTraj" + args[0] + ".path"

        PVKRB_traj_path = ExpName + "Data/PVKRBTraj" + args[0] + ".txt"
        PVKCP_traj_path = ExpName + "Data/PVKCPTraj" + args[0] + ".txt"
        PVKHJB_traj_path = ExpName + "Data/PVKHJBTraj" + args[0] + ".txt"
        ZSC_traj_path = ExpName + "Data/ZSCTraj" + args[0] + ".txt"
        OE_traj_path = ExpName + "Data/OETraj" + args[0] + ".txt"
        CP_traj_path = ExpName + "Data/CPTraj" + args[0] + ".txt"
        ZMP_traj_path = ExpName + "Data/ZMPTraj" + args[0] + ".txt"

    else:
        raise RuntimeError("Only one path file avaiable at one time!")

    # Load trajectories one by one
    robot_traj_tot = Traj_Loader(robot_traj_path)
    robot_traj = []
    time_count = 0;
    for robot_traj_i in robot_traj_tot:
        # Duration time
        delta_t_index = robot_traj_i.find('\t')
        if time_count ==0:
            delta_t = float(robot_traj_i[0:delta_t_index])
        time_count = 1
        robot_traj_i = robot_traj_i[delta_t_index + 1:]
        DOF_index = robot_traj_i.find('\t')
        DOF = int(robot_traj_i[0:DOF_index])
        robot_traj_i = robot_traj_i[DOF_index+1:]
        end_index = robot_traj_i.find('\n')
        robot_traj_i = robot_traj_i[0:end_index]
        robot_traj_i = robot_traj_i.split(" ")
        robot_traj_i = String_List_to_Number_List(robot_traj_i[0:DOF],"float")
        robot_traj.append(robot_traj_i)

    PVKRB_traj = String_List_to_Number_List(Traj_Loader(PVKRB_traj_path),"float")
    PVKCP_traj = String_List_to_Number_List(Traj_Loader(PVKCP_traj_path),"float")
    PVKHJB_traj = String_List_to_Number_List(Traj_Loader(PVKHJB_traj_path),"float")
    ZSC_traj = String_List_to_Number_List(Traj_Loader(ZSC_traj_path),"float")
    OE_traj = String_List_to_Number_List(Traj_Loader(OE_traj_path),"float")
    CP_traj = String_List_to_Number_List(Traj_Loader(CP_traj_path),"float")
    ZMP_traj = String_List_to_Number_List(Traj_Loader(ZMP_traj_path),"float")

    return delta_t, DOF, robot_traj, PVKRB_traj, PVKCP_traj, PVKHJB_traj, ZSC_traj, OE_traj, CP_traj, ZMP_traj

def Contact_Link_Reader(File_Name, Path_Name):
    # This function is used to read-in the formation of the certain link and its associated contact points
    # The output of this function is a dictionary of several keys with multiple list values
    # File_Name = "./User_File/Contact_Link.txt"
    Contact_Link_Dictionary = dict()
    # The format of this function should be an integet with a list of contact points
    Link_Number_i = -1
    File_Path_Name = Path_Name +File_Name
    with open(File_Path_Name) as Txt_File:
        Txt_File_Str = Txt_File.read().splitlines()
        Dictionary_Value_Add_Flag = 0
        for Txt_File_Str_i in Txt_File_Str:
            if "Link" in Txt_File_Str_i:
                # This indicates a contact link number
                # Then the next few points would be the local coordinates of the contact extremities
                Txt_File_Str_i = Txt_File_Str_i.translate(None, 'Link')     # This step is to get the link number out
                Link_Number_i = int(Txt_File_Str_i)
                Contact_Link_Dictionary[Link_Number_i] = []
                continue
            if Link_Number_i>=0:
                # Here the Txt_File_Str_i is a string of 'a, b, c' format
                Txt_File_Str_i = Txt_File_Str_i.split(",")
                del Txt_File_Str_i[-1]
                Txt_File_Flt_i = String_List_to_Number_List(Txt_File_Str_i,"float")
                Contact_Link_Dictionary[Link_Number_i].append(Txt_File_Flt_i)
    return Contact_Link_Dictionary

def Contact_Status_Reader(File_Name, Path_Name):
    # This function is used to read-in the formation of the certain link and its associated contact points
    # The output of this function is a dictionary of several keys with multiple list values
    # File_Name = "./User_File/Contact_Link.txt"
    Contact_Status_Dictionary = dict()
    # The format of this function should be an integet with a list of contact points
    Link_Number_i = -1
    File_Path_Name = Path_Name + File_Name
    with open(File_Path_Name) as Txt_File:
        Txt_File_Str = Txt_File.read().splitlines()
        Dictionary_Value_Add_Flag = 0
        for Txt_File_Str_i in Txt_File_Str:
            if "Link" in Txt_File_Str_i:
                # This indicates a contact link number
                # Then the next few points would be the local coordinates of the contact extremities
                Txt_File_Str_i = Txt_File_Str_i.translate(None, 'Link')     # This step is to get the link number out
                Link_Number_i = int(Txt_File_Str_i)
                Contact_Status_Dictionary[Link_Number_i] = []
                continue
            if Link_Number_i>=0:
                Contact_Status_Dictionary[Link_Number_i].append(int(Txt_File_Str_i))
    return Contact_Status_Dictionary

def Convex_Edge_Reader(File_Name, Path_Name):
    # This function is used to load in the Convex Edge for visualization
    # The format of this function should be an integer with a list of contact points
    File_Path_Name = Path_Name + File_Name
    Convex_Edge_List = []
    with open(File_Path_Name) as Txt_File:
        Txt_File_Str = Txt_File.read().splitlines()
        Dictionary_Value_Add_Flag = 0
        for Txt_File_Str_i in Txt_File_Str:
            Txt_File_Str_i = Txt_File_Str_i.split(" ")
            Txt_File_Flt_i = String_List_to_Number_List(Txt_File_Str_i,"float")
            Convex_Edge_List.append(Txt_File_Flt_i)
    return Convex_Edge_List

def Intersection_Reader(File_Name, Path_Name):
    File_Path_Name = Path_Name + File_Name
    Intersection_List = []
    with open(File_Path_Name) as Txt_File:
        Txt_File_Str = Txt_File.read().splitlines()
        Dictionary_Value_Add_Flag = 0
        for Txt_File_Str_i in Txt_File_Str:
            Txt_File_Str_i = Txt_File_Str_i.split(" ")
            Txt_File_Flt_i = String_List_to_Number_List(Txt_File_Str_i,"float")
            Intersection_List.append(Txt_File_Flt_i)
    return Intersection_List

def PIP_Traj_Reader(index):

    PIPList = []
    EdgeAList = []
    EdgeBList = []
    EdgeCOMList = []
    EdgexList = []
    EdgeyList = []
    EdgezList = []

    EdgeA_path = ExpName + "Data/EdgeATraj" + index + ".txt"
    with open(EdgeA_path,'r') as EdgeA_path_file:
        EdgeA_tot = EdgeA_path_file.readlines()
        for i in range(0, len(EdgeA_tot)):
            EdgeAList_i = []
            EdgeA_tot_i = EdgeA_tot[i].split(" ")
            EdgeAString = String_List_to_Number_List(EdgeA_tot_i[0:-1], "float")
            for j in range(0, len(EdgeAString)/3):
                EdgeAList_i.append(EdgeAString[3*j:3*j+3])
            EdgeAList.append(EdgeAList_i)

    EdgeB_path = ExpName + "Data/EdgeBTraj" + index + ".txt"
    with open(EdgeB_path,'r') as EdgeB_path_file:
        EdgeB_tot = EdgeB_path_file.readlines()
        for i in range(0, len(EdgeB_tot)):
            EdgeBList_i = []
            EdgeB_tot_i = EdgeB_tot[i].split(" ")
            EdgeBString = String_List_to_Number_List(EdgeB_tot_i[0:-1], "float")
            for j in range(0, len(EdgeBString)/3):
                EdgeBList_i.append(EdgeBString[3*j:3*j+3])
            EdgeBList.append(EdgeBList_i)

    EdgeCOM_path = ExpName + "Data/EdgeCOMTraj" + index + ".txt"
    with open(EdgeCOM_path,'r') as EdgeCOM_path_file:
        EdgeCOM_tot = EdgeCOM_path_file.readlines()
        for i in range(0, len(EdgeCOM_tot)):
            EdgeCOMList_i = []
            EdgeCOM_tot_i = EdgeCOM_tot[i].split(" ")
            EdgeCOMString = String_List_to_Number_List(EdgeCOM_tot_i[0:-1], "float")
            for j in range(0, len(EdgeCOMString)/3):
                EdgeCOMList_i.append(EdgeCOMString[3*j:3*j+3])
            EdgeCOMList.append(EdgeCOMList_i)

    Edgex_path = ExpName + "Data/EdgexTraj" + index + ".txt"
    with open(Edgex_path,'r') as Edgex_path_file:
        Edgex_tot = Edgex_path_file.readlines()
        for i in range(0, len(Edgex_tot)):
            EdgexList_i = []
            Edgex_tot_i = Edgex_tot[i].split(" ")
            EdgexString = String_List_to_Number_List(Edgex_tot_i[0:-1], "float")
            for j in range(0, len(EdgexString)/3):
                EdgexList_i.append(EdgexString[3*j:3*j+3])
            EdgexList.append(EdgexList_i)

    Edgey_path = ExpName + "Data/EdgeyTraj" + index + ".txt"
    with open(Edgey_path,'r') as Edgey_path_file:
        Edgey_tot = Edgey_path_file.readlines()
        for i in range(0, len(Edgey_tot)):
            EdgeyList_i = []
            Edgey_tot_i = Edgey_tot[i].split(" ")
            EdgeyString = String_List_to_Number_List(Edgey_tot_i[0:-1], "float")
            for j in range(0, len(EdgeyString)/3):
                EdgeyList_i.append(EdgeyString[3*j:3*j+3])
            EdgeyList.append(EdgeyList_i)

    Edgez_path = ExpName + "Data/EdgezTraj" + index + ".txt"
    with open(Edgez_path,'r') as Edgez_path_file:
        Edgez_tot = Edgez_path_file.readlines()
        for i in range(0, len(Edgez_tot)):
            EdgezList_i = []
            Edgez_tot_i = Edgez_tot[i].split(" ")
            EdgezString = String_List_to_Number_List(Edgez_tot_i[0:-1], "float")
            for j in range(0, len(EdgezString)/3):
                EdgezList_i.append(EdgezString[3*j:3*j+3])
            EdgezList.append(EdgezList_i)
    PIPList.append(EdgeAList)
    PIPList.append(EdgeBList)
    PIPList.append(EdgeCOMList)
    PIPList.append(EdgexList)
    PIPList.append(EdgeyList)
    PIPList.append(EdgezList)
    return PIPList

def PIP_Info_Reader(File_Name, Path_Name):
    # This function is used to load in the Convex Edge for visualization
    # The format of this function should be an integer with a list of contact points
    File_Path_Name = Path_Name + File_Name
    EdgeAList = []
    EdgeBList = []
    EdgeCOMList = []
    EdgexList = []
    EdgeyList = []
    EdgezList = []
    with open(File_Path_Name) as Txt_File:
        Txt_File_Str = Txt_File.read().splitlines()
        for i in range(0, len(Txt_File_Str)/6):
            EdgeAString = String_List_to_Number_List(Txt_File_Str[6*i].split(" "), "float")
            EdgeAList.append(EdgeAString)
            EdgeBString = String_List_to_Number_List(Txt_File_Str[6*i + 1].split(" "), "float")
            EdgeBList.append(EdgeBString)
            COMString = String_List_to_Number_List(Txt_File_Str[6*i + 2].split(" "), "float")
            EdgeCOMList.append(COMString)
            EdgexString = String_List_to_Number_List(Txt_File_Str[6*i + 3].split(" "), "float")
            EdgexList.append(EdgexString)
            EdgeyString = String_List_to_Number_List(Txt_File_Str[6*i + 4].split(" "), "float")
            EdgeyList.append(EdgeyString)
            EdgezString = String_List_to_Number_List(Txt_File_Str[6*i + 5].split(" "), "float")
            EdgezList.append(EdgezString)
    PIPList = []
    PIPList.append(EdgeAList)
    PIPList.append(EdgeBList)
    PIPList.append(EdgeCOMList)
    PIPList.append(EdgexList)
    PIPList.append(EdgeyList)
    PIPList.append(EdgezList)
    return PIPList

def Convex_Edges_Plot(sim_robot, convex_edges_list, vis):
    Convex_Edges_Number = len(convex_edges_list)/2
    COM_Pos = sim_robot.getCom()
    for i in range(0, Convex_Edges_Number):
        EdgeA = convex_edges_list[2*i]
        EdgeB = convex_edges_list[2*i + 1]
        # Three edges to be added: A->B, A -> COM, B-> COM
        Edge_Index = str(i)
        vis.add("Edge:" + Edge_Index, Trajectory([0, 1], [EdgeA, EdgeB]))
        vis.hideLabel("Edge:" + Edge_Index, True)
        vis.setAttribute("Edge:" + Edge_Index,'width', 5.0)

def Robot_Config_Plot(world, DOF, state_ref, contact_link_dictionary, convex_edges_list, delta_t=0.5):
    # This function is used to plot the robot motion
    # The optimized solution is used to plot the robot motion and the contact forces

    # Initialize the robot motion viewer
    robot_viewer = MyGLPlugin(world)
    # Here it is to unpack the robot optimized solution into a certain sets of the lists

    vis.pushPlugin(robot_viewer)
    vis.add("world", world)
    vis.show()

    sim_robot = world.robot(0)
    sim_robot.setConfig(state_ref[0:DOF])

    contact_link_list = contact_link_dictionary.keys()

    while vis.shown():
        # This is the main plot program
        vis.lock()
        sim_robot.setConfig(state_ref[0:DOF])
        # Convex_Edges_Plot(sim_robot, convex_edges_list, vis)
        # Robot_COM_Plot(sim_robot, vis)
        vis.unlock()
        time.sleep(delta_t)

def PIP_Subplot(i, EdgeA, EdgeB, EdgeCOM, Edgex, Edgey, Edgez, COM_Pos, vis):
    scale = 0.25

    Edge_Index = str(i)
    vis.add("PIPEdge:" + Edge_Index, Trajectory([0, 1], [EdgeA, EdgeB]))
    vis.hideLabel("PIPEdge:" + Edge_Index, True)
    vis.setAttribute("PIPEdge:" + Edge_Index,'width', 7.5)

    vis.add("PIPEdgeCOM:" + Edge_Index, Trajectory([0, 1], [COM_Pos, EdgeCOM]))
    vis.hideLabel("PIPEdgeCOM:" + Edge_Index, True)
    vis.setAttribute("PIPEdgeCOM:" + Edge_Index,'width', 7.5)
    vis.setColor("PIPEdgeCOM:" + Edge_Index, 65.0/255.0, 199.0/255.0, 244.0/255.0, 1.0)

    # Local Coordinates
    Edgex_i = [ 0.0, 0.0, 0.0]
    Edgex_i[0] = EdgeCOM[0] + scale * Edgex[0]
    Edgex_i[1] = EdgeCOM[1] + scale * Edgex[1]
    Edgex_i[2] = EdgeCOM[2] + scale * Edgex[2]
    vis.add("PIPEdgex:" + Edge_Index, Trajectory([0, 1], [EdgeCOM, Edgex_i]))
    vis.hideLabel("PIPEdgex:" + Edge_Index, True)
    vis.setAttribute("PIPEdgex:" + Edge_Index,'width', 7.5)
    vis.setColor("PIPEdgex:" + Edge_Index, 1.0, 0.0, 0.0, 1.0)

    Edgey_i = [ 0.0, 0.0, 0.0]
    Edgey_i[0] = EdgeCOM[0] + scale * Edgey[0]
    Edgey_i[1] = EdgeCOM[1] + scale * Edgey[1]
    Edgey_i[2] = EdgeCOM[2] + scale * Edgey[2]
    vis.add("PIPEdgey:" + Edge_Index, Trajectory([0, 1], [EdgeCOM, Edgey_i]))
    vis.hideLabel("PIPEdgey:" + Edge_Index, True)
    vis.setAttribute("PIPEdgey:" + Edge_Index,'width', 7.5)
    vis.setColor("PIPEdgey:" + Edge_Index, 155.0/255.0, 244.0/255.0, 66.0/255.0, 1.0)

    Edgez_i = [ 0.0, 0.0, 0.0]
    Edgez_i[0] = EdgeCOM[0] + scale * Edgez[0]
    Edgez_i[1] = EdgeCOM[1] + scale * Edgez[1]
    Edgez_i[2] = EdgeCOM[2] + scale * Edgez[2]
    # print Edgez_i
    vis.add("PIPEdgez:" + Edge_Index, Trajectory([0, 1], [EdgeCOM, Edgez_i]))
    vis.hideLabel("PIPEdgez:" + Edge_Index, True)
    vis.setAttribute("PIPEdgez:" + Edge_Index,'width', 7.5)
    vis.setColor("PIPEdgez:" + Edge_Index, 68.0/255.0, 65.0/255.0, 244.0/255.0, 1.0)

def PIP_Remove(i, vis):
    # This function is used to remove the previous PIP to make sure that no residual PIPs exist.
    Edge_Index = str(i)
    vis.remove("PIPEdge:" + Edge_Index)
    vis.remove("PIPEdgeCOM:" + Edge_Index)
    vis.remove("PIPEdgex:" + Edge_Index)
    vis.remove("PIPEdgey:" + Edge_Index)
    vis.remove("PIPEdgez:" + Edge_Index)

def PIP_Config_Plot(world, DOF, state_ref, pips_list, CPFlag, delta_t=0.5):
    # This function is used to plot the robot motion
    # The optimized solution is used to plot the robot motion and the contact forces

    # Initialize the robot motion viewer
    robot_viewer = MyGLPlugin(world)
    # Here it is to unpack the robot optimized solution into a certain sets of the lists

    vis.pushPlugin(robot_viewer)
    vis.add("world", world)
    vis.show()

    sim_robot = world.robot(0)
    sim_robot.setConfig(state_ref[0:DOF])

    InfeasiFlag = 0

    while vis.shown():
        # This is the main plot program
        vis.lock()
        sim_robot.setConfig(state_ref[0:DOF])
        PIPs_Number = len(pips_list[0])
        COM_Pos = sim_robot.getCom()
        EdgeAList = pips_list[0]
        EdgeBList = pips_list[1]
        EdgeCOMList = pips_list[2]
        EdgexList = pips_list[3]
        EdgeyList = pips_list[4]
        EdgezList = pips_list[5]

        for i in range(0, PIPs_Number):
            EdgeA = EdgeAList[i]
            EdgeB = EdgeBList[i]
            EdgeCOM = EdgeCOMList[i]
            Edgey = EdgexList[i]
            Edgez = EdgeyList[i]
            Edgex = EdgezList[i]
            PIP_Subplot(i, EdgeA, EdgeB, EdgeCOM, Edgex, Edgey, Edgez, COM_Pos, vis)

        Robot_COM_Plot(sim_robot, vis)

        if CPFlag is 1:
            try:
                h = ConvexHull(EdgeAList)
            except:
                InfeasiFlag = 1
        if InfeasiFlag is 0 and CPFlag is 1:
            h = ConvexHull(EdgeAList)
            hrender = draw_hull.PrettyHullRenderer(h)
            vis.add("blah", h)
            vis.setDrawFunc("blah", my_draw_hull)

        vis.unlock()
        time.sleep(delta_t)

def Robot_Traj_Plot(world, DOF, state_traj, contact_link_dictionary, delta_t=0.5):
    # Initialize the robot motion viewer
    robot_viewer = MyGLPlugin(world)
    # Here it is to unpack the robot optimized solution into a certain sets of the lists

    vis.pushPlugin(robot_viewer)
    vis.add("world", world)
    vis.show()

    sim_robot = world.robot(0)
    contact_link_list = contact_link_dictionary.keys()

    while vis.shown():
        # This is the main plot program
        for config_i in state_traj:
            vis.lock()
            sim_robot.setConfig(config_i[0:DOF])
            # Robot_COM_Plot(sim_robot, vis)
            vis.unlock()
            time.sleep(delta_t)

def Robot_COM_Plot(sim_robot, vis):
    COMPos_start = sim_robot.getCom()
    COMPos_end = COMPos_start[:]
    COMPos_end[2] = COMPos_end[2] - 7.50
    vis.add("COM", Trajectory([0, 1], [COMPos_start, COMPos_end]))
    vis.hideLabel("COM",True)
    vis.setColor("COM", 0.0, 204.0/255.0, 0.0, 1.0)
    vis.setAttribute("COM",'width', 7.5)
    # print COMPos_start

def GhostPlot(vis, NumberOfPoses, StepIndex, PoseSepIndex):
    FrameIndex = int(math.floor(StepIndex/PoseSepIndex))-1
    if FrameIndex == -1:
        return
    for i in range(0, FrameIndex):
        vis.hide("Ghost" + str(i),hidden=False)

def AssumedContacts(sim_robot, contact_link_dictionary, contact_status_dictionary):
    # This function is used to get all the active contacts out
    ActiveContacts = []
    for i in contact_link_dictionary.keys():
        for j in contact_status_dictionary[i]:
            if j is 1:
                ActiveContact = sim_robot.link(i).getWorldPosition(contact_link_dictionary[i][j])
                ActiveContacts.append(ActiveContact)
    return ActiveContacts

def FailureIndicator(Traj):
    FailureFlag = False
    FailureIndex = 0
    for TrajIndex in range(0, len(Traj)):
        if(Traj[TrajIndex]>0):
            FailureFlag = True
            FailureIndex = TrajIndex + 1
            return FailureFlag, FailureIndex
    return FailureFlag, FailureIndex

def Traj_Vis(world, DOF, robot_traj, PIP_traj, CPFlag, delta_t=0.5):
    # Initialize the robot motion viewer
    robot_viewer = MyGLPlugin(world)

    # Here it is to unpack the robot optimized solution into a certain sets of the lists
    vis.pushPlugin(robot_viewer)
    vis.add("world", world)
    vis.show()

    # ipdb.set_trace()
    state_traj = robot_traj[0]
    sim_robot = world.robot(0)
    EdgeAList = PIP_traj[0]
    EdgeBList = PIP_traj[1]
    EdgeCOMList = PIP_traj[2]
    EdgexList = PIP_traj[3]
    EdgeyList = PIP_traj[4]
    EdgezList = PIP_traj[5]

    TotalTrajLength = len(state_traj)
    EffectiveTrajLength = len(EdgeAList)
    RedundantTrajLength = TotalTrajLength - EffectiveTrajLength

    NumberOfConfigs = len(state_traj)

    TransScale = 0.35

    PVKRBcolor = [0, 0.4470, 0.7410]
    PVKCPcolor = [0.6588, 0.6588, 0.1961]
    PVKHJBcolor = [0.4940, 0.1840, 0.5560]
    ZSCcolor = [0.3010, 0.7450, 0.9330]
    OEcolor = [0.6350, 0.0780, 0.1840]
    CPcolor = [0.25, 0.25, 0.25]
    ZMPcolor = [0.75, 0, 0.75]

    if CPFlag is 2:
        state_traj = robot_traj[0]
        PVKRB_traj = robot_traj[1]
        PVKCP_traj = robot_traj[2]
        PVKHJB_traj = robot_traj[3]
        ZSC_traj = robot_traj[4]
        OE_traj = robot_traj[5]
        CP_traj = robot_traj[6]
        ZMP_traj = robot_traj[7]

        PVKRBFlag, PVKRBIndex = FailureIndicator(PVKRB_traj)
        PVKCPFlag, PVKCPIndex = FailureIndicator(PVKCP_traj)
        PVKHJBFlag, PVKHJBIndex = FailureIndicator(PVKHJB_traj)
        ZSCFlag, ZSCIndex = FailureIndicator(ZSC_traj)
        OEFlag, OEIndex = FailureIndicator(OE_traj)
        CPFlag, CPIndex = FailureIndicator(CP_traj)
        ZMPFlag, ZMPIndex = FailureIndicator(ZMP_traj)

        Flags = [PVKRBFlag, PVKCPFlag, PVKHJBFlag, ZSCFlag, OEFlag, CPFlag, ZMPFlag]
        Indices = [PVKRBIndex, PVKCPIndex, PVKHJBIndex, ZSCIndex, OEIndex, CPIndex, ZMPIndex]

        if PVKRBFlag is True:
            GhostPose = state_traj[PVKRBIndex + RedundantTrajLength]
            vis.add("PVKRBGhost", GhostPose)
            vis.hide("PVKRBGhost")
            vis.setColor("PVKRBGhost", PVKRBcolor[0], PVKRBcolor[1], PVKRBcolor[2], TransScale)

        if PVKCPFlag is True:
            GhostPose = state_traj[PVKCPIndex + RedundantTrajLength]
            vis.add("PVKCPGhost", GhostPose)
            vis.hide("PVKCPGhost")
            vis.setColor("PVKCPGhost", PVKCPcolor[0], PVKCPcolor[1], PVKCPcolor[2], TransScale)

        if PVKHJBFlag is True:
            GhostPose = state_traj[PVKHJBIndex + RedundantTrajLength]
            vis.add("PVKHJBGhost", GhostPose)
            vis.hide("PVKHJBGhost")
            vis.setColor("PVKHJBGhost", PVKHJBcolor[0], PVKHJBcolor[1], PVKHJBcolor[2], TransScale)

        if ZSCFlag is True:
            GhostPose = state_traj[ZSCIndex + RedundantTrajLength]
            vis.add("ZSCGhost", GhostPose)
            vis.hide("ZSCGhost")
            vis.setColor("ZSCGhost", ZSCcolor[0], ZSCcolor[1], ZSCcolor[2], TransScale)

        if OEFlag is True:
            GhostPose = state_traj[OEIndex + RedundantTrajLength]
            vis.add("OEGhost", GhostPose)
            vis.hide("OEGhost")
            vis.setColor("OEGhost", OEcolor[0], OEcolor[1], OEcolor[2], TransScale)

        if CPFlag is True:
            GhostPose = state_traj[CPIndex + RedundantTrajLength]
            vis.add("CPGhost", GhostPose)
            vis.hide("CPGhost")
            vis.setColor("CPGhost", CPcolor[0], CPcolor[1], CPcolor[2], TransScale)

        if ZMPFlag is True:
            GhostPose = state_traj[ZMPIndex + RedundantTrajLength]
            vis.add("ZMPGhost", GhostPose)
            vis.hide("ZMPGhost")
            vis.setColor("ZMPGhost", ZMPcolor[0], ZMPcolor[1], ZMPcolor[2], TransScale)

    while vis.shown():
        # This is the main plot program
        InfeasiFlag = 0
        for i in range(0, EffectiveTrajLength):
            vis.lock()
            config_i = state_traj[i + RedundantTrajLength]
            sim_robot.setConfig(config_i[0:DOF])
            COM_Pos = sim_robot.getCom()
            Robot_COM_Plot(sim_robot, vis)
            EdgeAList_i = EdgeAList[i]
            EdgeBList_i = EdgeBList[i]
            EdgeCOMList_i = EdgeCOMList[i]
            EdgexList_i = EdgexList[i]
            EdgeyList_i = EdgeyList[i]
            EdgezList_i = EdgezList[i]

            for j in range(0, len(EdgeAList_i)):
                EdgeA = EdgeAList_i[j]
                EdgeB = EdgeBList_i[j]
                EdgeCOM = EdgeCOMList_i[j]
                Edgex = EdgexList_i[j]
                Edgey = EdgeyList_i[j]
                Edgez = EdgezList_i[j]
                PIP_Subplot(j, EdgeA, EdgeB, EdgeCOM, Edgex, Edgey, Edgez, COM_Pos, vis)

            # # For the 7 fall indicators
            # if (PVKRBFlag is True) and (i>=PVKRBIndex):
            #     vis.hide("PVKRBGhost", False)
            # if (PVKCPFlag is True) and (i>=PVKCPIndex):
            #     vis.hide("PVKCPGhost", False)
            # if (PVKHJBFlag is True) and (i>=PVKHJBIndex):
            #     vis.hide("PVKHJBGhost", False)
            # if (ZSCFlag is True) and (i>=ZSCIndex):
            #     vis.hide("ZSCGhost", False)
            # if (OEFlag is True) and (i>OEIndex):
            #     vis.hide("OEGhost", False)
            # if (CPFlag is True) and (i>CPIndex):
            #     vis.hide("CPGhost", False)
            # if (ZMPFlag is True) and (i>ZMPIndex):
            #     vis.hide("ZMPGhost", False)
            # print Flags
            # print Indices

            if CPFlag is 1 or 2:
                try:
                    h = ConvexHull(EdgeAList_i)
                except:
                    InfeasiFlag = 1
                if InfeasiFlag is 0:
                    h = ConvexHull(EdgeAList_i)
                    hrender = draw_hull.PrettyHullRenderer(h)
                    vis.add("blah", h)
                    vis.setDrawFunc("blah", my_draw_hull)
                else:
                    print "Input Contact Polytope Infeasible!"
            vis.unlock()
            time.sleep(delta_t)

            for j in range(0, len(EdgeAList_i)):
                PIP_Remove(j, vis)

            if((CPFlag is 1 or 2) and (InfeasiFlag is 0)):
                vis.remove("blah")

        # if (PVKRBFlag is True):
        #     vis.hide("PVKRBGhost")
        # if (PVKCPFlag is True):
        #     vis.hide("PVKCPGhost")
        # if (PVKHJBFlag is True):
        #     vis.hide("PVKHJBGhost")
        # if (ZSCFlag is True):
        #     vis.hide("ZSCGhost")
        # if (OEFlag is True):
        #     vis.hide("OEGhost")
        # if (CPFlag is True):
        #     vis.hide("CPGhost")
        # if (ZMPFlag is True):
        #     vis.hide("ZMPGhost")
        print "End"

def ContactVerticesAppender(CP, point):
    # This function is used to append vertex to the end of contact vertices
    eps = 1e-8
    CPNo = len(CP);
    for j in range(0, len(point)):
        # Given point should be a 3-dimensional point
        point_j = point[j]
        ValueList = [0] * CPNo
        for i in range(0, CPNo):
            CP_i = CP[i]
            ValueList_i_x = CP_i[0] - point_j[0]
            ValueList_i_y = CP_i[1] - point_j[1]
            ValueList_i_z = CP_i[2] - point_j[2]
            ValueList_i = ValueList_i_x * ValueList_i_x + ValueList_i_y * ValueList_i_y + ValueList_i_z * ValueList_i_z
            ValueList[i] = ValueList_i
        if min(ValueList)>eps:
            CP.append(point_j)

def ContactVerticesGenerator(ConvexEdgesList):
    # This function is used to gather the contact vertices together.
    ConvexVertices = []
    ConvexVertices.append(ConvexEdgesList[0])
    for i in range(0, len(ConvexEdgesList)-1):
        ContactVerticesAppender(ConvexVertices, [ConvexEdgesList[i+1]])
    return ConvexVertices

def COM2IntersectionPlot(i, vis, COM_Pos, Intersection):
    # This function is used to plot PIP from COM to intersections
    Edge_Index = str(i)
    vis.add("PIPEdgefromCOM:" + Edge_Index, Trajectory([0, 1], [COM_Pos, Intersection]))
    vis.hideLabel("PIPEdgefromCOM:" + Edge_Index, True)
    vis.setAttribute("PIPEdgefromCOM:" + Edge_Index,'width', 7.5)
    vis.setColor("PIPEdgefromCOM:" + Edge_Index, 65.0/255.0, 199.0/255.0, 244.0/255.0, 1.0)

def main(*arg):

    Exp_Name = arg[0]
    VisFlag = arg[1]

    # The default robot to be loaded is the HRP2 robot.
    Robot_Option = "../user/hrp2/"
    world = WorldModel()                    # WorldModel is a pre-defined class

    ipdb.set_trace()

    if ".path" in Exp_Name:
        # Then this means that we are given the robot trajectories for experimentations.
        Real_Exp_Name = copy.deepcopy(Exp_Name)
        Realer_Exp_Name = ""
        for str_i in Real_Exp_Name:
            if str_i == ".":
                break;
            Realer_Exp_Name = Realer_Exp_Name + str_i;
        delta_t, DOF, robot_traj, PVKRB_traj, PVKCP_traj, PVKHJB_traj, ZSC_traj, OE_traj, CP_traj, ZMP_traj = Traj_Loader_fn(Realer_Exp_Name)
        # The next step is to load in robot's XML file
        XML_path = ExpName + "Envi/Envi" + Realer_Exp_Name + ".xml"
        result = world.readFile(XML_path)         # Here result is a boolean variable indicating the result of this loading operation
        if not result:
            raise RuntimeError("Unable to load model " + fn)
        Contact_Link_Dictionary = Contact_Link_Reader("ContactLink.txt", ExpName+"User/")
        Contact_Status_Dictionary = Contact_Status_Reader("InitContact.txt", ExpName+"User/")

        # There are five Vis flag options:
        """
            * VisFlag = 0 : Pure Animation
            * VisFlag = 1 : Pure Animation with ePIPs and edges
            * VisFlag = 2 : Pure Animation with eCH and ePIPs
            * VisFlag = 3 : Pure Animation With eCH, ePIPs, and Ghosts
        """

        if(VisFlag == 0):
            Robot_Traj_Plot(world, DOF, robot_traj, Contact_Link_Dictionary, delta_t)
        else:
            if(VisFlag == 1):
                PIPList = PIP_Traj_Reader(Realer_Exp_Name)
                CPFlag = 0
                Traj_Vis(world, DOF, robot_traj, PIPList, CPFlag, delta_t)
            else:
                if VisFlag is 2:
                    PIPList = PIP_Traj_Reader(Realer_Exp_Name)
                    CPFlag = 1
                    Traj_Vis(world, DOF, robot_traj, PIPList, CPFlag, delta_t)
                else:
                    if VisFlag is 3:
                        PIPList = PIP_Traj_Reader(Realer_Exp_Name)
                        CPFlag = 2
                        Traj_Vis(world, DOF, [robot_traj, PVKRB_traj, PVKCP_traj, PVKHJB_traj, ZSC_traj, OE_traj, CP_traj, ZMP_traj], PIPList, CPFlag, delta_t)
                    else:
                        raise RuntimeError("Flag value not feasible for trajectory visualization!")
    else:
        # In this case, what we have is a config
        # The following function can be used in two ways: the first way is to load the Config_Init.config file while the second way is to load two
        DOF, Config_Init, Velocity_Init = State_Loader_fn(Exp_Name + ".config", Robot_Option)
        # DOF, Config_Init, Velocity_Init = State_Loader_fn("Init_Config.txt", "Init_Velocity.txt", Robot_Option)
        # DOF, Config_Init, Velocity_Init = State_Loader_fn("Opt_Init_Config.txt", "Opt_Init_Velocity.txt", Robot_Option)

        # State_Writer_fn(Config_Init, "Init_Config_from_txt.config", Robot_Option)
        # State_Writer_fn(Config_Init, Velocity_Init, "Inn_Config.txt", "Inn_Velo.txt", Robot_Option)

        # According to the initial condition of the robot contact status, a basic optimization may have to be contacted to enforce the initial constraints.
        # Now it is the validation of the feasibility of the given initial condition
        State_Init = Config_Init + Velocity_Init
        Convex_Edges_List = Convex_Edge_Reader(Exp_Name + "CHEdges.txt", Robot_Option);

        delta_t = 0.5

        if(PIP_Flag == 0):
            Robot_Config_Plot(world, DOF, State_Init, Contact_Link_Dictionary, Convex_Edges_List)
        else:
            if PIP_Flag is 1:
                # In this case, we have to load in another sets of information for visualization
                PIPList = PIP_Info_Reader(Exp_Name + "PIPs.txt", Robot_Option)
                PIP_Config_Plot(world, DOF, State_Init, PIPList, 0)
            else:
                if PIP_Flag is 2:
                    PIPList = PIP_Info_Reader(Exp_Name + "PIPs.txt", Robot_Option)
                    PIP_Config_Plot(world, DOF, State_Init, PIPList, 1)
                else:
                    if PIP_Flag is 3:
                        # This function is specific for the four visualization of robot contact polytope
                        # Four plots should be shown alternatively

                        # 1. Contact Polytope with Edges
                        # 2. Contact Polytope with PIPsExpName
                        # 3. Contact Polytope with EPIPs
                        # 4. E Contact Polytope with EPIPs

                        Convex_Edges_List = Convex_Edge_Reader(Exp_Name + "CHEdges.txt", Robot_Option)
                        ContactVerticies = ContactVerticesGenerator(Convex_Edges_List)

                        Intersection_List = Intersection_Reader(Exp_Name + "Intersections.txt", Robot_Option)
                        PIPList = PIP_Info_Reader(Exp_Name + "PIPs.txt", Robot_Option)

                        # Initialize the robot motion viewer
                        robot_viewer = MyGLPlugin(world)
                        # Here it is to unpack the robot optimized solution into a certain sets of the lists

                        vis.pushPlugin(robot_viewer)
                        vis.add("world", world)
                        vis.show()

                        sim_robot = world.robot(0)
                        sim_robot.setConfig(State_Init[0:DOF])

                        COM_Pos = sim_robot.getCom()
                        # t_orbital = ...
                        # q_orbital = ...
                        # vis.add("q_orbital",q_orbital)
                        # vis.hide("q_orbital")
                        # vis.setColor("q_orbital",1,0,0,0.5)

                        while vis.shown():
                            # This is the main plot program
                            vis.lock()
                            # sim_robot.setConfig(Config_Init[0:DOF])

                            # if vis.animationTime() < t_orbital:
                            #     vis.hide("q_orbital")
                            # else:
                            #     vis.hide("q_orbital",False)

                            if robot_viewer.mode_no == 1:
                                # Here we plot the full contact polytope
                                h = ConvexHull(ContactVerticies)
                                hrender = draw_hull.PrettyHullRenderer(h)
                                vis.add("blah", h)
                                vis.setDrawFunc("blah", my_draw_hull)

                            if robot_viewer.mode_no == 2:
                                h = ConvexHull(ContactVerticies)
                                hrender = draw_hull.PrettyHullRenderer(h)
                                vis.add("blah", h)
                                vis.setDrawFunc("blah", my_draw_hull)

                                # Also edges
                                Convex_Edges_Plot(sim_robot, Convex_Edges_List, vis)

                                # Then intersection from COM
                                for i in range(0, len(Intersection_List)):
                                    COM2IntersectionPlot(i, vis, COM_Pos, Intersection_List[i])
                            if robot_viewer.mode_no == 3:
                                h = ConvexHull(ContactVerticies)
                                hrender = draw_hull.PrettyHullRenderer(h)
                                vis.add("blah", h)
                                vis.setDrawFunc("blah", my_draw_hull)

                                # Also edges
                                Convex_Edges_Plot(sim_robot, Convex_Edges_List, vis)

                                PIPs_Number = len(PIPList[0])
                                EdgeAList = PIPList[0]
                                EdgeBList = PIPList[1]
                                EdgeCOMList = PIPList[2]
                                EdgexList = PIPList[3]
                                EdgeyList = PIPList[4]
                                EdgezList = PIPList[5]

                                for i in range(0, PIPs_Number):
                                    EdgeA = EdgeAList[i]
                                    EdgeB = EdgeBList[i]
                                    EdgeCOM = EdgeCOMList[i]
                                    Edgey = EdgexList[i]
                                    Edgez = EdgeyList[i]
                                    Edgex = EdgezList[i]
                                    PIP_Subplot(i, EdgeA, EdgeB, EdgeCOM, Edgex, Edgey, Edgez, COM_Pos, vis)

                            if robot_viewer.mode_no == 4:
                                # Effective Polytope with Effective PIPs
                                PIPs_Number = len(PIPList[0])

                                EdgeAList = PIPList[0]
                                # EdgeBList = PIPList[1]
                                # for i in EdgeBList:
                                #     EdgeAList.append(i);
                                h = ConvexHull(EdgeAList)
                                hrender = draw_hull.PrettyHullRenderer(h)
                                vis.add("blah", h)
                                vis.setDrawFunc("blah", my_draw_hull)

                                EdgeAList = PIPList[0]
                                EdgeBList = PIPList[1]
                                EdgeCOMList = PIPList[2]
                                EdgexList = PIPList[3]
                                EdgeyList = PIPList[4]
                                EdgezList = PIPList[5]

                                for i in range(0, PIPs_Number):
                                    EdgeA = EdgeAList[i]
                                    EdgeB = EdgeBList[i]
                                    EdgeCOM = EdgeCOMList[i]
                                    Edgey = EdgexList[i]
                                    Edgez = EdgeyList[i]
                                    Edgex = EdgezList[i]
                                    PIP_Subplot(i, EdgeA, EdgeB, EdgeCOM, Edgex, Edgey, Edgez, COM_Pos, vis)
                            vis.unlock()
                            print "mode number: " + str(robot_viewer.mode_no)
                            time.sleep(delta_t)

                            # Robot_COM_Plot(sim_robot, vis)

                            # Now we should delete the imposed object
                            if robot_viewer.mode_no == 1:
                                vis.remove("blah")
                            if robot_viewer.mode_no == 2:
                                vis.remove("blah")
                                EdgeNo = len(Convex_Edges_List)/2
                                for i in range(0, EdgeNo):
                                    Edge_Index = str(i)
                                    vis.remove("Edge:" + Edge_Index)
                                for i in range(0, len(Intersection_List)):
                                    Edge_Index = str(i)
                                    vis.remove("PIPEdgefromCOM:" + Edge_Index)
                            if robot_viewer.mode_no == 3:
                                vis.remove("blah")
                                EdgeNo = len(Convex_Edges_List)/2
                                for i in range(0, EdgeNo):
                                    Edge_Index = str(i)
                                    vis.remove("Edge:" + Edge_Index)
                                EdgeAList = PIPList[0]
                                for i in range(0, len(EdgeAList)):
                                    Edge_Index = str(i)
                                    vis.remove("PIPEExpNamedge:" + Edge_Index)
                                    vis.remove("PIPEdgeCOM:" + Edge_Index)
                                    vis.remove("PIPEdgex:" + Edge_Index)
                                    vis.remove("PIPEdgey:" + Edge_Index)
                                    vis.remove("PIPEdgez:" + Edge_Index)
                            if robot_viewer.mode_no == 4:
                                vis.remove("blah")
                                EdgeAList = PIPList[0]
                                for i in range(0, len(EdgeAList)):
                                    Edge_Index = str(i)
                                    vis.remove("PIPEdge:" + Edge_Index)
                                    vis.remove("PIPEdgeCOM:" + Edge_Index)
                                    vis.remove("PIPEdgex:" + Edge_Index)
                                    vis.remove("PIPEdgey:" + Edge_Index)
                                    vis.remove("PIPEdgez:" + Edge_Index)

                            robot_viewer.mode_no = robot_viewer.mode_no + 1
                            if robot_viewer.mode_no == 5:
                                robot_viewer.mode_no = 1
                            robot_viewer.mode_no = 3

if __name__ == "__main__":
    # main("InitConfig", 2)
    main("3.path", 3)
