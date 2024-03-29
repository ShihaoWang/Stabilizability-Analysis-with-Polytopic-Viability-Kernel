#include <ctime>
#include "CommonHeader.h"
#include "NonlinearOptimizerInfo.h"
#include <omp.h>
#include "Control/PathController.h"
#include "Simulation/WorldSimulation.h"
#include <ode/ode.h>

static bool InnerSimulation(const std::string & FolderPath, ViabilityKernelInfo & VKObj)
{
  RobotWorld world;
  SimGUIBackend Backend(&world);
  WorldSimulation& Sim = Backend.sim;

  /* 0. Load the XML World file */
  string XMLFileStr = FolderPath + "/Sample.xml";
  const char* XMLFile = XMLFileStr.c_str();    // Here we must give abstract path to the file
  if(!Backend.LoadAndInitSim(XMLFile))
  {
    std::cerr<<"Sample XML file path does not exist!"<<endl;
    return true;
  }

  /* 1. Load the Contact Link file */
  const std::string UserFilePath = FolderPath + "/user/hrp2/";
  const std::string ContactLinkPath = UserFilePath + "ContactLink.txt";
  int NumberOfContactPoints;
  std::vector<LinkInfo> RobotLinkInfo = ContactInfoLoader(ContactLinkPath, NumberOfContactPoints);

  /* 2. Load the Contact Status file */
  const std::string ContactStatusPath = UserFilePath + "InitContact.txt";
  std::vector<ContactStatusInfo> RobotContactInfo = ContactStatusInfoLoader(ContactStatusPath);

  /* 3. Generation of robot's configuration and World XML file */
  Robot SimRobot = *world.robots[0];
  std::vector<double> InitRobotConfig(SimRobot.q.size()), InitRobotVelocity(SimRobot.q.size()), ZeroRobotVelocity(SimRobot.q.size());
  SignedDistanceFieldInfo SDFInfo = InitEnviGenerator(SimRobot, RobotLinkInfo, RobotContactInfo, UserFilePath, InitRobotConfig);
  RobotConfigWriter(InitRobotConfig, UserFilePath, "InitConfig.config");

  // Here we would like to give it one last chance, with which the disturbance is given with impulse.
  Config RobotConfigNew(InitRobotConfig);
  SimRobot.UpdateConfig(RobotConfigNew);
  SimRobot.UpdateGeometry();
  SimRobot.dq = ZeroRobotVelocity;

  //  Given the optimized result to be the initial state
  Config InitRobotConfigNew(InitRobotConfig);
  Sim.world->robots[0]->UpdateConfig(InitRobotConfigNew);
  Sim.world->robots[0]->dq = ZeroRobotVelocity;

  Config InitRobotVelocityNew(ZeroRobotVelocity);
  Sim.controlSimulators[0].oderobot->SetConfig(InitRobotConfigNew);
  Sim.controlSimulators[0].oderobot->SetVelocities(InitRobotVelocityNew);

  /* 4. Projected Inverted Pendulum Plot*/
  std::vector<Vector3> ActContactPositions, ActVelocities;
  std::vector<Matrix> ActJacobians;
  std::vector<double> ActDists;
  std::vector<int> ActStatus = ActContactNJacobian(SimRobot, RobotLinkInfo, RobotContactInfo, ActContactPositions, ActVelocities, ActDists, ActJacobians, SDFInfo);

  std::vector<Vector3> CPVertex, CPEdgeA, CPEdgeB;
  std::vector<FacetInfo> FacetInfoObj = ContactHullGeneration(ActContactPositions, CPVertex, CPEdgeA, CPEdgeB);      // This function output is only used for visualization purpose.
  ConvexEdgesWriter(FacetInfoObj, UserFilePath, "InitConfigCHEdges.txt");

  Vector3 InitCOM = SimRobot.GetCOM();
  int FailureFlag;
  std::vector<PIPInfo> PIPTotal = ContactEdgesGenerationSP(ActContactPositions, ActVelocities, ActDists, ActStatus, InitCOM, InitCOM, FailureFlag);   // SP denotes SP projection approach.
  PIPsWriter(PIPTotal, UserFilePath, "InitConfigPIPs.txt");
  std::vector<Vector3> FullPIPInters = FullPIPInterCal(FacetInfoObj, InitCOM);
  IntersectionsWriter(FullPIPInters, UserFilePath, "InitConfigIntersections.txt");

  // This part is to simulate robot's world file.
  RobotWorld worldAct;
  SimGUIBackend BackendAct(&worldAct);
  WorldSimulation& SimAct = BackendAct.sim;

  int FileIndex = FileIndexFinder();
  XMLFileStr = FolderPath + "/build/Envi" + std::to_string(FileIndex) + ".xml";
  XMLFile = XMLFileStr.c_str();    // Here we must give abstract path to the file
  if(!BackendAct.LoadAndInitSim(XMLFile))
  {
    std::cerr<<"Envi" + std::to_string(FileIndex)<<".xml file path does not exist!"<<endl;
    return true;
  }

  double t_imp = 1.0;
  // Here first we would like to run the simulation for a certain amount of time to ensure the contact is well established.
  double  dt          = 0.025;
  string stateTrajFile = "stateTraj" + std::to_string(FileIndex) + ".path";
  const char *stateTrajFile_Name = stateTrajFile.c_str();
  while(SimAct.time <= t_imp)
  {
    SimAct.Advance(dt);
    SimAct.UpdateModel();
    BackendAct.DoStateLogging_LinearPath(0, stateTrajFile_Name);
  }

  /* 5. Internal Experimentation */
  SimulationTest(SimAct, VKObj, RobotLinkInfo, RobotContactInfo, SDFInfo, BackendAct, dt, FileIndex);
  return true;
}

int main()
{
  string ViabilityKernelPath = "/home/motion/Desktop/VKACC/build/2.25/";
  bool VKFastFlag = false;
  ViabilityKernelInfo VKObj = ViabilityKernelDataLoader(ViabilityKernelPath, VKFastFlag);
  for (int i = 0; i < 1001; i++)
  {
    if(i%5==0)
    {
      string mv_command = "mv -f *Traj*.* ./Data";
      const char *mv_command_str = mv_command.c_str();
      std::system(mv_command_str);

      mv_command = "mv -f Specs*.* ./Data";
      mv_command_str = mv_command.c_str();
      std::system(mv_command_str);

      mv_command = "mv -f Envi*.* ./Envi";
      mv_command_str = mv_command.c_str();
      std::system(mv_command_str);
    }
    // Simulation loop at each time
    std::string FolderPath = "/home/motion/Desktop/Stabilizability-Analysis-with-Polytopic-Viability-Kernel";
    bool InitializationFlag = false;
    while (InitializationFlag == false)
    {
      InitializationFlag = InnerSimulation(FolderPath, VKObj);
    }
  }
  return 0;
}
