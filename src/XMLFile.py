# This function is used to generate robot world file for simulation experimentation.
import sys, os, time
import ipdb
import copy
import math
import numpy as np
from klampt.model import ik,coordinates,config,trajectory,collide
import klampt

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

def Configuration_Writer_fn(DOF, Config):
    ConfigString = "\"" + str(DOF) + "\t"
    for i in Config:
        ConfigString+=str(i) + " "
    ConfigString+="\"/>"
    return ConfigString

def RotXYAngle(Axis):
    # This function is used to calculate Euler Angles from rotation axis
    Axis_x = Axis[0]
    Axis_y = Axis[1]
    Axis_z = Axis[2]
    RotY = math.asin(Axis_x)
    RotX = -math.atan2(Axis_y, Axis_z)
    return RotX, RotY

def main():
    Robot_Option = "../user/hrp2/"
    AxesPath = Robot_Option + "Axes.txt"
    Axes = []
    LinkIndices = []
    RotX = []
    RotY = []
    with open(AxesPath,'r') as Axes_file:
        AxesTot = Axes_file.readlines()
        AxesNo = len(AxesTot)/4
        for i in range(0, AxesNo):
            AxesInit = 4 * i
            AxesEnd = AxesInit + 4
            AxesInfo = AxesTot[AxesInit:AxesEnd]
            Axis_i = []
            for i in range(0,3):
                Axis_val = AxesInfo[i]
                Axis_val = float(Axis_val[:-2])
                Axis_i.append(Axis_val)
            LinkIndex = AxesInfo[3]
            LinkIndex = LinkIndex[:-1]
            LinkIndex = int(LinkIndex)
            RotX_i, RotY_i = RotXYAngle(Axis_i)
            Axes.append(Axis_i)
            LinkIndices.append(LinkIndex)
            RotX.append(RotX_i)
            RotY.append(RotY_i)
    CenterPoint = []
    CenterPointPath = Robot_Option + "AvgContacts.txt"
    with open(CenterPointPath,'r') as CenterPoint_file:
        CenterPointTotal = CenterPoint_file.readlines()
        CenterPointNo = len(CenterPointTotal)/3
        for i in range(0, CenterPointNo):
            CPInit = 3 * i
            CPEnd = CPInit + 3
            CPInfo = CenterPointTotal[CPInit:CPEnd]
            CP_i = []
            for i in range(0,3):
                CP_val = CPInfo[i]
                CP_val = float(CP_val[:-2])
                CP_i.append(CP_val)
            CenterPoint.append(CP_i)

    # YawPath = Robot_Option + "YawAngles.txt"
    # Yaw = []
    # with open(YawPath,'r') as Yaw_file:
    #     Yaws = Yaw_file.readlines()
    #     for Yaw_i in Yaws:
    #         Yaw_i = Yaw_i.replace("\n", "")
    #         Yaw.append(Yaw_i)

    # alright now the job is to write them into a world file
    ConfigName = Robot_Option + "Config4XML.config"
    DOF, Config = Configuration_Loader_fn(ConfigName)

    XMLHeader = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<world>\n"
    ConfigInfo = Configuration_Writer_fn(DOF, Config)
    RobotHeader = "  <robot name=\"HRP2\" file=\"/home/motion/Klampt-examples/data/robots/hrp2.rob\" config=" + ConfigInfo
    TerrainHeader =""
    for i in range(0, len(CenterPoint)):
        TerrainHeader+="\n  <terrain file=\"/home/motion/Klampt-examples/data/terrains/block.off\" scale=\"0.1 0.1 0.35\" translation=\""
        TerrainHeader+= str(CenterPoint[i][0]) + " " + str(CenterPoint[i][1]) + " " + str(CenterPoint[i][2]) + "\""
        TerrainHeader+=" rotateX= \"" + str(RotX[i]) + " \" rotateY=\"" + str(RotY[i]) + "\">\n     <display color=\"0.4 0.3 0.2\"/>\n  </terrain>"
    SimulationHeasder = "\n  <simulation>\n    <globals adaptiveTimeStepping=\"1\" />"
    for i in range(0, len(CenterPoint)):
        SimulationHeasder+="\n    <robot index=\"0\" body=\"" + str(LinkIndices[i]) + "\">"
        SimulationHeasder+="\n      <geometry kFriction=\"1.0\" kRestitution=\"0.0\" padding=\"0.005\" stiffness=\"80000\" damping=\"20000\" />\n    </robot>"
    SimulationHeasder+="\n  </simulation> \n</world>\n"
    TotalFile = XMLHeader + RobotHeader + TerrainHeader + SimulationHeasder

    ExpIndexPath = "../build/FileIndex.txt"
    ExpIndexNumber = 0
    with open(ExpIndexPath,'r') as ExpIndexFile:
        ExpIndex = ExpIndexFile.readlines()
        ExpIndex = ExpIndex[0].replace("\n", "")
        ExpIndexNumber = int(ExpIndex)

    XMLFileName = "../build/" + "Envi" + str(ExpIndexNumber) + ".xml"
    with open(XMLFileName, 'w') as file:
        file.write(TotalFile)
    # ipdb.set_trace()
    world = klampt.WorldModel()
    res = world.readFile(XMLFileName)
    if not res:
        raise RuntimeError("Unable to load model "+fn)
    collider = collide.WorldCollider(world)
    TerrainNo = len(collider.terrains)
    ActLinkList = []
    for i in range(0, TerrainNo):
        Iterator = collider.robotTerrainCollisions(world.robot(0), world.terrain(i))
        for CollisionPair_i in Iterator:
            ActLinkList.append(CollisionPair_i[0].getIndex())
    Status = Contact_Status_Reader("InitContact.txt", Robot_Option)
    ActLinks = Status.keys()
    TerrainCollisionFlag = 0
    for i in ActLinkList:
        if i not in ActLinks:
            TerrainCollisionFlag = 1
    # The last step is to write TerrainCollisionFlag
    TerrainCollisionName = Robot_Option + "TerrainCollision.txt"
    with open(TerrainCollisionName, 'w') as file:
        file.write(str(TerrainCollisionFlag))
if __name__ == "__main__":
    main()
