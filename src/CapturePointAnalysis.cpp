#include "RobotInfo.h"
#include "CommonHeader.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <iterator>

// This function is used to calculate objective trajectories for PVK-CP and CP.
static std::vector<double> linspace(const double & start_in, const double & end_in, const int & num_in)
{
  std::vector<double> linspaced;
  double start = static_cast<double>(start_in);
  double end = static_cast<double>(end_in);
  double num = static_cast<double>(num_in);
  if (num == 0) { return linspaced; }
  if (num == 1)
    {
      linspaced.push_back(start);
      return linspaced;
    }
  double delta = (end - start) / (num - 1);
  for(int i=0; i < num-1; ++i)
    {
      linspaced.push_back(start + delta * i);
    }
  linspaced.push_back(end);
  return linspaced;
}

static std::vector<double> getVertexIndices(std::string const& pointLine)
{
  std::istringstream iss(pointLine);

  return std::vector<double>{
    std::istream_iterator<double>(iss),
    std::istream_iterator<double>()
  };
}

// static void qTrajNqdotTrajLoader(const string & UserPath, const int& FileIndex, std::vector<Config> & qTraj, std::vector<Config> & qdotTraj, std::vector<double> & COMVelx, std::vector<double> & COMVely, std::vector<double> & COMVelz)
static void qTrajNqdotTrajLoader(const string & UserPath, const int& FileIndex, std::vector<Config> & qTraj, std::vector<Config> & qdotTraj)
{
  const string qTrajFile = "qTrajAct" + std::to_string(FileIndex) + ".txt";
  const char *qTrajFile_Name = qTrajFile.c_str();
  const string qdotTrajFile = "qdotTrajAct" + std::to_string(FileIndex) + ".txt";
  const char *qdotTrajFile_Name = qdotTrajFile.c_str();

  string str_line;
  string ConfigFilePath = UserPath + qTrajFile_Name;
  ifstream ConfigInfofile (ConfigFilePath);
  if (ConfigInfofile.is_open())
  {
    while (getline(ConfigInfofile, str_line) )
    {
      std::vector<double> RobotConfig = getVertexIndices(str_line);
      Config RobotConfig_i(RobotConfig);
      qTraj.push_back(RobotConfig_i);
    }
    ConfigInfofile.close();
  }
  else
  {
    std::cerr<<"Wrong! State config txt cannot be found!"<<endl;
  }

  string VelocityFilePath = UserPath + qdotTrajFile_Name;
  ifstream VelocityInfofile (VelocityFilePath);
  if (VelocityInfofile.is_open())
  {
    while (getline (VelocityInfofile, str_line) )
    {
      // std::cout << str_line << "\n";
      std::vector<double> RobotVelocity = getVertexIndices(str_line);
      // VectorPrintResult(RobotVelocity);
      Config RobotVelocity_i(RobotVelocity);
      qdotTraj.push_back(RobotVelocity_i);
    }
    VelocityInfofile.close();
  }
  else
  {
    std::cerr<<"Wrong! State velocity txt cannot be found!"<<endl;
  }

  return;
}

static void ObjTrajWriter(const std::vector<double> & Traj, const int & ExpIndex, const string & FallDetector)
{
  // This function is used to generate the ROC curve
  string ObjTrajFile = FallDetector + "Traj" + std::to_string(ExpIndex) + ".txt";
  const char *ObjTrajFile_Name = ObjTrajFile.c_str();
  std::ofstream ObjWriter;
  ObjWriter.open(ObjTrajFile_Name);
  for (int i = 0; i < Traj.size(); i++)
  {
    ObjWriter<<std::to_string(Traj[i])<<"\n";
  }
  ObjWriter.close();
}

static void CPEvaluation(const int & FileIndex, Robot& SimRobot, ViabilityKernelInfo& VKObj, std::vector<LinkInfo> & RobotLinkInfo, std::vector<ContactStatusInfo> & RobotContactInfo, SignedDistanceFieldInfo & SDFInfo, const std::vector<Config> & qTraj, const std::vector<Config> & qdotTraj)
{
  // This function takes in the given trajectory of q and qdot and output for each method.
  double dt = 0.025;        // This should be set in advance according to simulation.

  int EdgeNumber = 4;
  double mu = 1.0;

  // Seven objective trajectories
  const int StepNo = qTraj.size();
  std::vector<double> PVKRBTraj(StepNo),        PVKCPTraj(StepNo),        PVKHJBTraj(StepNo);
  std::vector<double> ZSCTraj(StepNo),          OETraj(StepNo),           CPTraj(StepNo),           ZMPTraj(StepNo);

  // Centroidal Informatioin only contains robot's centroidal position and velocity.
  std::vector<double> COMx(StepNo),             COMy(StepNo),             COMz(StepNo);
  std::vector<double> COMVelx(StepNo),          COMVely(StepNo),          COMVelz(StepNo);

  std::vector<const char*> EdgeFileNames;
  string fEdgeAFile = "EdgeATraj" + std::to_string(FileIndex) + ".txt";                 const char *fEdgeAFile_Name = fEdgeAFile.c_str();
  string fEdgeBFile = "EdgeBTraj" + std::to_string(FileIndex) + ".txt";                 const char *fEdgeBFile_Name = fEdgeBFile.c_str();
  string fEdgeCOMFile = "EdgeCOMTraj" + std::to_string(FileIndex) + ".txt";             const char *fEdgeCOMFile_Name = fEdgeCOMFile.c_str();
  string fEdgexTrajFile = "EdgexTraj" + std::to_string(FileIndex) + ".txt";             const char *fEdgexTrajFile_Name = fEdgexTrajFile.c_str();
  string fEdgeyTrajFile = "EdgeyTraj" + std::to_string(FileIndex) + ".txt";             const char *fEdgeyTrajFile_Name = fEdgeyTrajFile.c_str();
  string fEdgezTrajFile = "EdgezTraj" + std::to_string(FileIndex) + ".txt";             const char *fEdgezTrajFile_Name = fEdgezTrajFile.c_str();

  EdgeFileNames.push_back(fEdgeAFile_Name);         EdgeFileNames.push_back(fEdgeBFile_Name);         EdgeFileNames.push_back(fEdgeCOMFile_Name);
  EdgeFileNames.push_back(fEdgexTrajFile_Name);     EdgeFileNames.push_back(fEdgeyTrajFile_Name);     EdgeFileNames.push_back(fEdgezTrajFile_Name);

  int StepIndex = 0;
  for (int i = 0; i < qTraj.size(); i++)        // For the sake of Centroidal Acceleration for ZMP
  {
    std::printf("Index %d\n", StepIndex);

    if (StepIndex ==46)
    {
      int a = 1;
    }

    Config qNow = qTraj[i];
    Config qdotNow = qdotTraj[i];
    SimRobot.UpdateConfig(qNow);
    SimRobot.dq = qdotNow;

    /* Robot's COMPos and COMVel */
    Vector3 COMPos(0.0, 0.0, 0.0), COMVel(0.0, 0.0, 0.0), COMAcc(0.0, 0.0, 0.0);
    CentroidalState(SimRobot, COMPos, COMVel);
    COMx[StepIndex] = COMPos.x;             COMy[StepIndex] = COMPos.y;               COMz[StepIndex] = COMPos.z;
    COMVelx[StepIndex] = COMVel.x;          COMVely[StepIndex] = COMVel.y;            COMVelz[StepIndex] = COMVel.z;

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

    // Capture Point which is a 2D versino of PVK-CP
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

    StepIndex = StepIndex + 1;
  }

  // // Here the job is to write down the centroidal trajectories.
  // ObjTrajWriter(COMx, FileIndex, "COMx");
  // ObjTrajWriter(COMy, FileIndex, "COMy");
  // ObjTrajWriter(COMz, FileIndex, "COMz");
  // ObjTrajWriter(COMVelx, FileIndex, "COMVelx");
  // ObjTrajWriter(COMVely, FileIndex, "COMVely");
  // ObjTrajWriter(COMVelz, FileIndex, "COMVelz");
  //
  // ObjTrajWriter(PVKRBTraj, FileIndex, "PVKRB");
  // ObjTrajWriter(PVKCPTraj, FileIndex, "PVKCP");
  // ObjTrajWriter(PVKHJBTraj, FileIndex, "PVKHJB");
  //
  // ObjTrajWriter(ZSCTraj, FileIndex, "ZSC");
  //
  // ObjTrajWriter(OETraj, FileIndex, "OE");
  // ObjTrajWriter(CPTraj, FileIndex, "CP");
  // ObjTrajWriter(ZMPTraj, FileIndex, "ZMP");

  return;
}

void CapturePointAnalysis(Robot & SimRobot, ViabilityKernelInfo & VKObj, std::vector<LinkInfo> & RobotLinkInfo, std::vector<ContactStatusInfo> & RobotContactInfo, SignedDistanceFieldInfo & SDFInfo)
{
  // This function is use to generate data analysis for experimentation trajectories.
  string UserPath = "/home/motion/Desktop/Stabilizability-Analysis-with-Polytopic-Viability-Kernel/build/Data/";
  for (int i = 483; i < 666; i++)
  {
    int FileIndex = i + 1;
    std::vector<Config> qTraj, qdotTraj;
    qTrajNqdotTrajLoader(UserPath, FileIndex, qTraj, qdotTraj);
    CPEvaluation(FileIndex, SimRobot, VKObj, RobotLinkInfo, RobotContactInfo, SDFInfo, qTraj, qdotTraj);
  }
  return;
}
