// This function is used to extensively simulate the result from the four failure metric
#include "RobotInfo.h"
#include "CommonHeader.h"
#include <ode/ode.h>
#include "Control/PathController.h"
#include "Control/JointTrackingController.h"

void SimulationTest(WorldSimulation & Sim, ViabilityKernelInfo& VKObj, std::vector<LinkInfo> & RobotLinkInfo, std::vector<ContactStatusInfo> & RobotContactInfo, SignedDistanceFieldInfo & SDFInfo, SimGUIBackend & Backend, const double & dt, const int & FileIndex)
{
  /* Simulation parameters */
  int     EdgeNumber      = 4;
  double  mu              = 1.0;
  double  t_last          = 2.0 * dt;
  double  t_cur           = Sim.time;
  double  t_impul         = t_cur + t_last;                      // The impulse lasts for 0.5s.
  double  t_final         = 4.0;                                  // The simulation lasts for 3.5s.
  int     StepNo          = round(t_final/dt);
  t_final+=Sim.time;

  std::random_device rd;
  std::mt19937 gen(rd());

  double ImpFx = 3500.0;
  double ImpFy = 3500.0;
  double ImpFz = 1000.0;

  std::uniform_real_distribution<> ImpXdis(-ImpFx, ImpFx);
  std::uniform_real_distribution<> ImpYdis(-ImpFy, ImpFy);
  std::uniform_real_distribution<> ImpZdis(-ImpFz, ImpFz);

  double Fx_t = ImpXdis(gen);
  double Fy_t = ImpYdis(gen);
  double Fz_t = ImpZdis(gen);
  Vector3 F_t(Fx_t, Fy_t, Fz_t);

  /* Override the default controller with a PolynomialPathController */
  auto NewControllerPtr = std::make_shared<PolynomialPathController>(*Sim.world->robots[0]);
  Sim.SetController(0, NewControllerPtr);
  NewControllerPtr->SetConstant(Sim.world->robots[0]->q);

  std::vector<Vector3> ContactPositionRef;
  // This function is used to get the robot's current active end effector position and Jacobian matrices.
  for (int i = 0; i < RobotLinkInfo.size(); i++)
  {
    for (int j = 0; j < RobotLinkInfo[i].LocalContacts.size(); j++)
    {
      switch (RobotContactInfo[i].LocalContactStatus[j])
      {
        case 1:
        {
          // This means that current contact is active and we should keep its location and Jacobian.
          Vector3 LinkiPjPos;
          Sim.world->robots[0]->GetWorldPosition(RobotLinkInfo[i].LocalContacts[j], RobotLinkInfo[i].LinkIndex, LinkiPjPos);
          ContactPositionRef.push_back(LinkiPjPos);
        }
        break;
        default:
        break;
      }
    }
  }

  /* Simulation Trajectory */
  // These three are used to save the trajectory of the desired robot's properties.
  std::vector<Config> qDesTraj,       qdotDesTraj;
  qDesTraj.reserve(StepNo);           qdotDesTraj.reserve(StepNo);
  qDesTraj.push_back(Sim.world->robots[0]->q);
  qdotDesTraj.push_back(Sim.world->robots[0]->dq);

  // These two save the trajectory of robot's actual properties.
  std::vector<Config> qActTraj,       qdotActTraj;
  qActTraj.reserve(StepNo);           qdotActTraj.reserve(StepNo);
  qActTraj.push_back(Sim.world->robots[0]->q);
  qdotActTraj.push_back(Sim.world->robots[0]->dq);

  // Centroidal Informatioin only contains robot's centroidal position and velocity.
  std::vector<double> COMx(StepNo),             COMy(StepNo),             COMz(StepNo);
  std::vector<double> COMVelx(StepNo),          COMVely(StepNo),          COMVelz(StepNo);

  std::vector<double> KETraj(StepNo);

  // Seven objective trajectories
  std::vector<double> PVKRBTraj(StepNo),        PVKCPTraj(StepNo),        PVKHJBTraj(StepNo);
  std::vector<double> ZSCTraj(StepNo),          OETraj(StepNo),           CPTraj(StepNo),           ZMPTraj(StepNo);

  int NumberOfActEndEffectorInit;               // This variable describes the number of active end effectors!
  int NumberOfContactPoints;                    // This variable describes the number of total contact points!
  ContactNumberFinder(RobotContactInfo, NumberOfActEndEffectorInit, NumberOfContactPoints);
  SpecsWriter(*Sim.world->robots[0], t_final, t_last, F_t, NumberOfActEndEffectorInit, FileIndex);

  int DOF = Sim.world->robots[0]->q.size();

  std::vector<const char*> EdgeFileNames, CentroidalFileNames, ObjectiveNames, StateTrajNames;
  string fEdgeAFile = "EdgeATraj" + std::to_string(FileIndex) + ".txt";                 const char *fEdgeAFile_Name = fEdgeAFile.c_str();
  string fEdgeBFile = "EdgeBTraj" + std::to_string(FileIndex) + ".txt";                 const char *fEdgeBFile_Name = fEdgeBFile.c_str();
  string fEdgeCOMFile = "EdgeCOMTraj" + std::to_string(FileIndex) + ".txt";             const char *fEdgeCOMFile_Name = fEdgeCOMFile.c_str();
  string fEdgexTrajFile = "EdgexTraj" + std::to_string(FileIndex) + ".txt";             const char *fEdgexTrajFile_Name = fEdgexTrajFile.c_str();
  string fEdgeyTrajFile = "EdgeyTraj" + std::to_string(FileIndex) + ".txt";             const char *fEdgeyTrajFile_Name = fEdgeyTrajFile.c_str();
  string fEdgezTrajFile = "EdgezTraj" + std::to_string(FileIndex) + ".txt";             const char *fEdgezTrajFile_Name = fEdgezTrajFile.c_str();

  EdgeFileNames.push_back(fEdgeAFile_Name);         EdgeFileNames.push_back(fEdgeBFile_Name);         EdgeFileNames.push_back(fEdgeCOMFile_Name);
  EdgeFileNames.push_back(fEdgexTrajFile_Name);     EdgeFileNames.push_back(fEdgeyTrajFile_Name);     EdgeFileNames.push_back(fEdgezTrajFile_Name);

  const string COMxFile = "COMxTraj" + std::to_string(FileIndex) + ".txt";              const char *COMxFile_Name = COMxFile.c_str();
  const string COMyFile = "COMyTraj" + std::to_string(FileIndex) + ".txt";              const char *COMyFile_Name = COMyFile.c_str();
  const string COMzFile = "COMzTraj" + std::to_string(FileIndex) + ".txt";              const char *COMzFile_Name = COMzFile.c_str();
  const string COMVelxFile = "COMVelxTraj" + std::to_string(FileIndex) + ".txt";        const char *COMVelxFile_Name = COMVelxFile.c_str();
  const string COMVelyFile = "COMVelyTraj" + std::to_string(FileIndex) + ".txt";        const char *COMVelyFile_Name = COMVelyFile.c_str();
  const string COMVelzFile = "COMVelzTraj" + std::to_string(FileIndex) + ".txt";        const char *COMVelzFile_Name = COMVelzFile.c_str();
  const string KEFile = "KETraj" + std::to_string(FileIndex) + ".txt";                  const char *KEFile_Name = KEFile.c_str();

  CentroidalFileNames.push_back(COMxFile_Name);     CentroidalFileNames.push_back(COMyFile_Name);     CentroidalFileNames.push_back(COMzFile_Name);
  CentroidalFileNames.push_back(COMVelxFile_Name);  CentroidalFileNames.push_back(COMVelyFile_Name);  CentroidalFileNames.push_back(COMVelzFile_Name);
  CentroidalFileNames.push_back(KEFile_Name);

  const string PVKRBTrajFile = "PVKRBTraj" + std::to_string(FileIndex) + ".txt";        const char *PVKRBTrajFile_Name = PVKRBTrajFile.c_str();
  const string PVKHJBTrajFile = "PVKHJBTraj" + std::to_string(FileIndex) + ".txt";      const char *PVKHJBTrajFile_Name = PVKHJBTrajFile.c_str();
  const string PVKCPTrajFile = "PVKCPTraj" + std::to_string(FileIndex) + ".txt";        const char *PVKCPTrajFile_Name = PVKCPTrajFile.c_str();
  const string ZSCTrajFile = "ZSCTraj" + std::to_string(FileIndex) + ".txt";            const char *ZSCTrajFile_Name = ZSCTrajFile.c_str();
  const string OETrajFile = "OETraj" + std::to_string(FileIndex) + ".txt";              const char *OETrajFile_Name = OETrajFile.c_str();
  const string CPTrajFile = "CPTraj" + std::to_string(FileIndex) + ".txt";              const char *CPTrajFile_Name = CPTrajFile.c_str();
  const string ZMPTrajFile = "ZMPTraj" + std::to_string(FileIndex) + ".txt";            const char *ZMPTrajFile_Name = ZMPTrajFile.c_str();

  ObjectiveNames.push_back(PVKRBTrajFile_Name);
  ObjectiveNames.push_back(PVKCPTrajFile_Name);
  ObjectiveNames.push_back(PVKHJBTrajFile_Name);
  ObjectiveNames.push_back(ZSCTrajFile_Name);
  ObjectiveNames.push_back(OETrajFile_Name);
  ObjectiveNames.push_back(CPTrajFile_Name);
  ObjectiveNames.push_back(ZMPTrajFile_Name);

  const string qActTrajFile = "qActTraj" + std::to_string(FileIndex) + ".txt";          const char *qActTrajFile_Name = qActTrajFile.c_str();
  const string qdotActTrajFile = "qdotActTraj" + std::to_string(FileIndex) + ".txt";    const char *qdotActTrajFile_Name = qdotActTrajFile.c_str();
  StateTrajNames.push_back(qActTrajFile_Name);
  StateTrajNames.push_back(qdotActTrajFile_Name);

  string stateTrajFile = "stateTraj" + std::to_string(FileIndex) + ".path";             const char *stateTrajFile_Name = stateTrajFile.c_str();

  /*
    Here we have three types of controller:
        1. Rigid-body controller
        2. Whole-body QP stabilizing controller
  */

  int ControllerType = 2;
  int StepIndex = 0;

  std::vector<double> qDes = Sim.world->robots[0]->q;   // This is commanded robot configuration to the controller.
  int QPStatus = 0;

  while(Sim.time < t_final)
  {
    switch (ControllerType)
    {
      case 1:
      {

      }
      break;
      default:
      {
        // A weird problem with actuators at ankles have been observed. Here we assume that these two motors behave in an ideal way.
        std::vector<double> RealVelocities(DOF);
        for (int i = 0; i < DOF; i++)
        {
          RealVelocities[i] = Sim.world->robots[0]->dq[i];
        }
        // Four actuators are compensated with ideal values.
        RealVelocities[10] = qdotDesTraj[qdotDesTraj.size()-1][10];
        RealVelocities[11] = qdotDesTraj[qdotDesTraj.size()-1][11];
        RealVelocities[16] = qdotDesTraj[qdotDesTraj.size()-1][16];
        RealVelocities[17] = qdotDesTraj[qdotDesTraj.size()-1][17];

        Sim.world->robots[0]->dq = RealVelocities;
        Config RealVelocitiesSet(RealVelocities);
        Sim.controlSimulators[0].oderobot->SetVelocities(RealVelocitiesSet);
      }
      break;
    }

    if(Sim.time<=t_impul)
    {
      // The impulse is given to the robot's torso
      dBodyAddForceAtPos(Sim.odesim.robot(0)->body(19), Fx_t, Fy_t, Fz_t, 0.0, 0.0, 0.0);
    }

    Robot SimRobot = *Sim.world->robots[0];

    /* Robot's COMPos and COMVel */
    Vector3 COMPos(0.0, 0.0, 0.0), COMVel(0.0, 0.0, 0.0), COMAcc(0.0, 0.0, 0.0);
    CentroidalState(SimRobot, COMPos, COMVel);
    COMx[StepIndex] = COMPos.x;             COMy[StepIndex] = COMPos.y;               COMz[StepIndex] = COMPos.z;
    COMVelx[StepIndex] = COMVel.x;          COMVely[StepIndex] = COMVel.y;            COMVelz[StepIndex] = COMVel.z;

    // Get robot's current kinetic energy
    KETraj[StepIndex] = SimRobot.GetKineticEnergy();

    std::vector<Vector3>  ActContactPositions, ActVelocities;        // A vector of Vector3 points
    std::vector<Matrix>   ActJacobians;       // A vector of Jacobian matrices
    std::vector<double>   ActDists;
    std::vector<int> ActStatus = ActContactNJacobian(SimRobot, RobotLinkInfo, RobotContactInfo, ActContactPositions, ActVelocities, ActDists, ActJacobians, SDFInfo);

    std::vector<Vector3> ConeShiftedUnits, ConeUnits;
    ConeUnitGenerator(ActContactPositions, SDFInfo, ConeShiftedUnits, ConeUnits, EdgeNumber, mu);

    /* 3. Failure Metric using PVK-HJB assumption*/
    double HJBObjective;
    std::vector<PIPInfo> PIPTotal = PIPGenerator(ActContactPositions, ActVelocities, ActDists, ActStatus, COMPos, COMVel, EdgeFileNames, VKObj, HJBObjective, dt);
    PVKHJBTraj[StepIndex] = HJBObjective;

    /* 1. Failure Metric using PVK-RB assumption*/
    double RBObjective = RBGenerator(PIPTotal);
    PVKRBTraj[StepIndex] = RBObjective;

    /* 2. Failure Metric using PVK-CP assumption*/
    double CPCEObjective = CPCEGenerator(PIPTotal);
    PVKCPTraj[StepIndex] = CPCEObjective;

    /* 4. Zero-Step-Capturability with Convex Optimization */
    double ZSCObjective = ZeroStepCapturabilityGenerator(ActContactPositions, ConeShiftedUnits, EdgeNumber, COMPos, COMVel);
    ZSCTraj[StepIndex] = ZSCObjective;

    // Projection to ground
    std::vector<Vector3> ProjActContactPositions, ZeroContactVelocities;
    ProjActContactPositions.reserve(ActContactPositions.size());
    ZeroContactVelocities.reserve(ActVelocities.size());
    double LowestHeight = 100.0;
    for (int j = 0; j < ActContactPositions.size(); j++)
    {
      if(LowestHeight>ActContactPositions[j].z)
      {
        LowestHeight = ActContactPositions[j].z;
      }
    }

    for (int j = 0; j < ActContactPositions.size(); j++)
    {
      Vector3 ProjActContact(ActContactPositions[j].x, ActContactPositions[j].y, LowestHeight);
      ProjActContactPositions.push_back(ProjActContact);

      Vector3 ZeroContactVelocity(0.0, 0.0, 0.0);
      ZeroContactVelocities.push_back(ZeroContactVelocity);
    }

    std::vector<double> PIPObj;
    double PVKHJBMargin = 0.0;
    double HJBSPObjective = 0.0;
    std::vector<PIPInfo> PIPSPTotal = PIPGeneratorAnalysis(ProjActContactPositions, ZeroContactVelocities, ActDists, ActStatus, COMPos, COMVel, VKObj, PIPObj, HJBSPObjective, PVKHJBMargin, dt);

    // Orbital Energy which is a 2D version of PVK-RB
    double OEMargin = 0.0;
    double OEObjective = RBGeneratorAnalysis(PIPSPTotal, OEMargin);
    OETraj[StepIndex] = OEObjective;

    // Capture Point which is a 2D version of PVK-CP
    double CPMargin = 0.0;
    double CPObjective = CPCEGeneratorAnalysis(PIPSPTotal, CPMargin);
    CPTraj[StepIndex] = CPObjective;

    switch (StepIndex)
    {
      case 0:
      {

      }
      break;
      default:
      {
        COMAcc.x = (COMVelx[StepIndex] - COMVelx[StepIndex-1])/dt;
        COMAcc.y = (COMVely[StepIndex] - COMVely[StepIndex-1])/dt;
        COMAcc.z = (COMVelz[StepIndex] - COMVelz[StepIndex-1])/dt;
      }
      break;
    }

    // ZMP
    double ZMPMargin = 0.0;
    double ZMPObjective = ZMPGeneratorAnalysis(PIPSPTotal, COMPos, COMAcc, ZMPMargin);
    ZMPTraj[StepIndex] = ZMPObjective;

    std::vector<double> FailureMetricVec = { RBObjective, CPCEObjective, HJBObjective, ZSCObjective, OEObjective, CPObjective, ZMPObjective};
    CentroidalFailureMetricWriter(COMPos, COMVel, KETraj[StepIndex], FailureMetricVec, CentroidalFileNames, ObjectiveNames);

    /*  Controller Input  */
    switch (ControllerType)
    {
      case 1:
      {
        // In this case, the robot's controller holds a constant initial configuration.
        std::printf("Using Controller 1: Rigid-body Controller!\n");

        qActTraj.push_back(SimRobot.q);
        qdotActTraj.push_back(SimRobot.dq);
      }
      break;
      case 2:
      {
        // In this case, the robot's controller would like to stabilize the robot with a QP controller.
        std::printf("Using Controller 2: QP Stabilizing Controller!\n");
        std::vector<double> qNew = StabilizingControllerContact(SimRobot, ActJacobians, ConeUnits, EdgeNumber, DOF, dt, qDesTraj, qdotDesTraj, qActTraj, qdotActTraj, RobotLinkInfo, RobotContactInfo, ContactPositionRef, ActContactPositions, ActVelocities, NumberOfContactPoints, StepIndex);
        qDes = qNew;
      }
      break;
      default:
      {
      }
      break;
    }

    TrajAppender(StateTrajNames[0], SimRobot.q, DOF);
    TrajAppender(StateTrajNames[1], SimRobot.dq, DOF);

    // Send the control command!
    Config qDesired(qDes);
    NewControllerPtr->SetConstant(qDesired);

    Sim.Advance(dt);
    Sim.UpdateModel();
    Backend.DoStateLogging_LinearPath(0, stateTrajFile_Name);
    StepIndex = StepIndex + 1;
  }
  return;
}
