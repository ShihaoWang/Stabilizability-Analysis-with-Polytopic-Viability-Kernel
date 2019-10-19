#include "NonlinearOptimizerInfo.h"
#include "CommonHeader.h"
#include "Planning/ConstraintChecker.h"
#include <bits/stdc++.h>
#include <string>     // std::string, std::to_string
#include <KrisLibrary/geometry/Conversions.h>
#include <KrisLibrary/geometry/MultiVolumeGrid.h>
#include <KrisLibrary/geometry/CollisionMesh.h>
#include <KrisLibrary/geometry/AnyGeometry.h>
#include <KrisLibrary/meshing/VolumeGrid.h>

// This function is used to generate randomized configurations whether certain constraints will be satisfied.

// Firstly,   the job is to randomly generate a configuration.
// Secondly,  the job is to optimize the configuration such that center of mass lies within support polygon.
// Thirdly,   the job is to try the self-collision checking to get rid of the collision configuration.
// Fourly,    the job is to relocate robot's center of mass position such that it is near origin while the minimum distance of its body part are at 0-height plane.

static Robot SimRobotObj;
static std::vector<LinkInfo> RobotLinkInfo;
static std::vector<ContactStatusInfo> RobotContactInfo;

static double PI = 3.1415926535897932;
static double ConeAngle = PI/6.0;        // Pi/6

static std::vector<double> ConfigSampler()
{
  // This function is used to sample the robot's configuration
  int DOF = SimRobotObj.q.size();
  std::vector<double> InitConfig(DOF);
  for (int i = 0; i < DOF; i++)
  {
    InitConfig[i] = 0.0;
  }
  // The first three coordinates are not important since eventually we would still have to relocate them.
  std::random_device rd;
  std::mt19937 gen(rd());
  for (int i = 3; i < DOF; i++)
  {
    // Configuration
    double qLow, qUpp;
    switch (i)
    {
      case 3:
      {
        // Yaw
        qLow = -PI;
        qUpp = PI;
      }
      break;
      case 4:
      {
        // Pitch
        qLow = -ConeAngle;
        qUpp = ConeAngle;
      }
      break;
      case 5:
      {
        // Roll
        qLow = -ConeAngle;
        qUpp = ConeAngle;
      }
      break;
      case 34:
      {
        // Left Hand
        qLow = -ConeAngle/2.0;
        qUpp = ConeAngle/2.0;
      }
      break;
      case 27:
      {
        // Right Hand
        qLow = -ConeAngle/2.0;
        qUpp = ConeAngle/2.0;
      }
      break;
      default:
      {
        qLow = SimRobotObj.qMin(i);
        qUpp = SimRobotObj.qMax(i);
      }
      break;
    }

    std::uniform_real_distribution<> qDis(qLow, qUpp);
    double q_i = qDis(gen);
    InitConfig[i] = q_i;
  }
  return InitConfig;
}

struct COMSPOpt: public NonlinearOptimizerInfo
{
  // This is the class struct for the optimizaiton of configuraiton variables

  COMSPOpt():NonlinearOptimizerInfo(){};

  // This struct inherits the NonlinearOptimizerInfo struct and we just need to defined the Constraint function
  static void ObjNConstraint(int    *Status, int *n,    double x[],
    int    *needF,  int *neF,  double F[],
    int    *needG,  int *neG,  double G[],
    char      *cu,  int *lencu,
    int    iu[],    int *leniu,
    double ru[],    int *lenru)
    {
      std::vector<double> x_vec(*n);
      for (int i = 0; i < *n; i++)
      {
        x_vec[i] = x[i];
      }
      std::vector<double> F_val = COMSPObjNCons(*n, *neF, x_vec);
      for (int i = 0; i < *neF; i++)
      {
        F[i] = F_val[i];
      }
    }
  void Solve(std::vector<double> &RobotConfig)
  {
    int StartType = 0;
    NonlinearProb.solve(StartType, neF, n, ObjAdd, ObjRow, ObjNConstraint,
      xlow, xupp, Flow, Fupp,
      x, xstate, xmul, F, Fstate, Fmul,
      nS, nInf, sumInf);
      for (int i = 0; i < n; i++)
      {
        RobotConfig[i] = x[i];
      }
      delete []x;      delete []xlow;   delete []xupp;
      delete []xmul;   delete []xstate;

      delete []F;      delete []Flow;   delete []Fupp;
      delete []Fmul;   delete []Fstate;
  }
  static std::vector<double> COMSPObjNCons(const int & nVar, const int & nObjNCons, const std::vector<double> & ConfigOpt)
  {
    // This funciton provides the constraint for the configuration variable

    std::vector<double> F(nObjNCons);
    Config ConfigOptNew(ConfigOpt);
    SimRobotObj.UpdateConfig(ConfigOptNew);     // Here both the SimRobot.q and robot frames have already been updated.

    double ConfigVia = 0.0;
    std::vector<double> ActiveContact;
    for (int i = 0; i < RobotLinkInfo.size(); i++)
    {
      int LinkiPNo = RobotLinkInfo[i].LocalContacts.size();
      int LinkiStatus = 0;
      for (int j = 0; j < LinkiPNo; j++)
      {
        switch (RobotContactInfo[i].LocalContactStatus[j])
        {
          case 0:
          break;
          case 1:
          {
            Vector3 LinkiPjPos;
            SimRobotObj.GetWorldPosition(RobotLinkInfo[i].LocalContacts[j], RobotLinkInfo[i].LinkIndex, LinkiPjPos);
            ActiveContact.push_back(LinkiPjPos.z);
            LinkiStatus = 1;
          }
          break;
          default:
          break;
        }
      }
      switch (LinkiStatus)
      {
        case 1:
        {
          int Link_index = RobotLinkInfo[i].LinkIndex;
          RobotLink3D Link_i = SimRobotObj.links[Link_index];
          double norm_x = Link_i.T_World.R.data[2][0];
          double norm_y = Link_i.T_World.R.data[2][1];
          double norm_z = Link_i.T_World.R.data[2][2];

          double TileAngle = PI/3.0;;
          switch (i)
          {
            case 0:
            {
              TileAngle = ConeAngle;
            }
            break;
            case 1:
            {
              TileAngle = ConeAngle;
            }
            break;
            default:
            {

            }
            break;
          }
          F[i+1] = norm_z - cos(TileAngle);
        }
        break;
        default:
        {
          F[i+1] = 0.0;
        }
        break;
      }
    }

    double ActiveLowestZ = *std::min_element(ActiveContact.begin(), ActiveContact.end());
    F[0] = ActiveLowestZ * ActiveLowestZ;

    // The objective is to let Center of Mass remains to be within SP as inner as possible.
    std::vector<Vector3> SPVertices;
    for (int i = 0; i < RobotLinkInfo.size(); i++)
    {
      int LinkiPNo = RobotLinkInfo[i].LocalContacts.size();
      for (int j = 0; j < LinkiPNo; j++)
      {
        switch (RobotContactInfo[i].LocalContactStatus[j])
        {
          case 0:
          break;
          case 1:
          {
            Vector3 LinkiPjPos;
            SimRobotObj.GetWorldPosition(RobotLinkInfo[i].LocalContacts[j], RobotLinkInfo[i].LinkIndex, LinkiPjPos);
            LinkiPjPos.z = 0.0;
            SPVertices.push_back(LinkiPjPos);
          }
          break;
          default:
          break;
        }
      }
    }
    Vector3 COM_Pos = SimRobotObj.GetCOM();
    int FacetFlag = 0;
    FacetInfo SPObj = FlatContactHullGeneration(SPVertices, FacetFlag);    // This is the support polygon
    COM_Pos.z = 0.0;
    F[5] = SPObj.ProjPoint2EdgeDist(COM_Pos) - 0.15;

    return F;
  }
};

static void CoMSPOpt(std::vector<double> & RobotConfig)
{
  // This function is used to optimize the robot configuraiton variable such that certain contact is active.
  COMSPOpt COMSPOptProblem;

  // Static Variable Substitution
  int n = SimRobotObj.q.size();
  int neF = 1;                  // Cost function how close projection COM to SP
  neF = neF + 5;                // Foot's normal axis direction lie within a cone.
  COMSPOptProblem.InnerVariableInitialize(n, neF);

  /*
    Initialize the bounds of variables
  */
  std::vector<double> xlow_vec(n), xupp_vec(n);

  for (int i = 0; i < n; i++)
  {
    // Configuration
    xlow_vec[i] = SimRobotObj.qMin(i);
    xupp_vec[i] = SimRobotObj.qMax(i);
  }
  COMSPOptProblem.VariableBoundsUpdate(xlow_vec, xupp_vec);

  /*
    Initialize the bounds of variables
  */
  std::vector<double> Flow_vec(neF), Fupp_vec(neF);
  Flow_vec[0] = 0;
  Fupp_vec[0] = 1e20;           // The default idea is that the objective should always be nonnegative

  Flow_vec[1] = 0;
  Fupp_vec[1] = 1e20;

  Flow_vec[2] = 0;
  Fupp_vec[2] = 1e20;

  Flow_vec[3] = 0;
  Fupp_vec[3] = 1e20;

  Flow_vec[4] = 0;
  Fupp_vec[4] = 1e20;

  Flow_vec[5] = 0;
  Fupp_vec[5] = 1e20;

  COMSPOptProblem.ConstraintBoundsUpdate(Flow_vec, Fupp_vec);

  /*
    Initialize the seed guess
  */
  COMSPOptProblem.SeedGuessUpdate(RobotConfig);

  /*
  Given a name of this problem for the output
  */
  COMSPOptProblem.ProblemNameUpdate("COMSPOptProblem", 0);

  COMSPOptProblem.NonlinearProb.setIntParameter("Iterations limit", 1000000);
  COMSPOptProblem.NonlinearProb.setIntParameter("Major iterations limit", 250);
  COMSPOptProblem.NonlinearProb.setIntParameter("Major print level", 0);
  COMSPOptProblem.NonlinearProb.setIntParameter("Minor print level", 0);
  /*
  ProblemOptions seting
  */
  // Solve with Finite-Difference
  COMSPOptProblem.ProblemOptionsUpdate(0, 3);
  COMSPOptProblem.Solve(RobotConfig);
  return;
}

std::vector<double> InitEnviGenerator(Robot & _SimRobotObj, const std::vector<LinkInfo> & _RobotLinkInfo, const std::vector<ContactStatusInfo> &  _RobotContactInfo, const string & UserFilePath, int & coLFlag)
{
  SimRobotObj = _SimRobotObj;
  RobotLinkInfo = _RobotLinkInfo;
  RobotContactInfo = _RobotContactInfo;

  std::vector<double> InitConfig;
  int SimuPassFlag = 0;
  while(SimuPassFlag ==0)
  {
    // The simulation flag has to be passed in order to proceed forward.
    bool TerrainCollisionFlag = true;
    while(TerrainCollisionFlag == true)
    {
      bool CollisionFlag = true;
      while(CollisionFlag == true)
      {
        coLFlag = 0;
        InitConfig = ConfigSampler();
        CoMSPOpt(InitConfig);
        Config RobotConfigNew(InitConfig);
        SimRobotObj.UpdateConfig(RobotConfigNew);
        SimRobotObj.UpdateGeometry();
        CollisionFlag = SimRobotObj.SelfCollision();
      }
      // If the robot's configuration can make it here, the next step is to construct a environment where robot's active contact
      RobotConfigWriter(InitConfig, UserFilePath, "Config4XML.config");

      // Get the normal vector directions out.
      std::vector<double> Axes;
      std::vector<Vector3> LocalContacts;
      std::vector<int> LinkIndices;
      std::vector<double> YawAngles;
      for (int i = 0; i < RobotLinkInfo.size(); i++)
      {
        switch (RobotContactInfo[i].LocalContactStatus[0])
        {
          case 0:
          break;
          case 1:
          {
            Vector3 LinkiPjPos;
            SimRobotObj.GetWorldPosition(RobotLinkInfo[i].AvgLocalContact, RobotLinkInfo[i].LinkIndex, LinkiPjPos);
            LocalContacts.push_back(LinkiPjPos);

            int Link_index = RobotLinkInfo[i].LinkIndex;
            LinkIndices.push_back(Link_index);

            RobotLink3D Link_i = SimRobotObj.links[Link_index];
            double norm_x = Link_i.T_World.R.data[2][0];
            double norm_y = Link_i.T_World.R.data[2][1];
            double norm_z = Link_i.T_World.R.data[2][2];
            Axes.push_back(norm_x);
            Axes.push_back(norm_y);
            Axes.push_back(norm_z);

            // One more thing needs to be done is to calcualte the Yaw angle for block rotation.
            Vector3 LinkiPjAPos, LinkiPjBPos;
            SimRobotObj.GetWorldPosition(RobotLinkInfo[i].LocalContacts[0], RobotLinkInfo[i].LinkIndex, LinkiPjAPos);
            SimRobotObj.GetWorldPosition(RobotLinkInfo[i].LocalContacts[1], RobotLinkInfo[i].LinkIndex, LinkiPjBPos);

            // Then we calculate the Yaw Angle
            double Dist_x = LinkiPjAPos.x - LinkiPjBPos.x;
            double Dist_y = LinkiPjAPos.y - LinkiPjBPos.y;
            double Yaw_i = atan2(Dist_y, Dist_x);
            YawAngles.push_back(Yaw_i);
          }
          break;
          default:
          break;
        }
      }
      // Now the job is to construct a RobotWorld XML file. The easy way is to use a Python file.
      RobotAxesWriter(Axes, LinkIndices, UserFilePath);
      // YawAngleWriter(YawAngles, UserFilePath);
      AvgContactWriter(LocalContacts, UserFilePath);
      std::system("cd ../src && python XMLFile.py");

      // Then we should read-in file to see whether a terrain collision happens.

      string str_line;
      string TerrainCollisionPath = UserFilePath + "TerrainCollision.txt";
      ifstream TerrainInfofile (TerrainCollisionPath);
      int TerrainCollisionFlag = 0;
      if (TerrainInfofile.is_open())
      {
        while (getline (TerrainInfofile, str_line) )
        {
          TerrainCollisionFlag = stoi(str_line);
        }
        TerrainInfofile.close();
      }
      else
      {
        std::cerr<<"Wrong! TerrainCollision.txt cannot be found!"<<endl;
      }

      switch (TerrainCollisionFlag)
      {
        case 0:
        {
          TerrainCollisionFlag = false;
        }
        break;
        default:
        {
          TerrainCollisionFlag = true;
        }
        break;
      }
    }

    int aaa = 1;

    // Now it is time to enable another world file for simulation test and collision test with environment objects.
    RobotWorld world;
    SimGUIBackend Backend(&world);
    WorldSimulation& Sim = Backend.sim;

    string FileIndexName = "FileIndex.txt";         // This file should be located in the "build" folder.
    ifstream FileIndexReader(FileIndexName);
    int FileIndex;
    string str_line;
    if (FileIndexReader.is_open())
    {
      while (getline (FileIndexReader,str_line) )
      {
        FileIndex = stoi(str_line);
      }
      FileIndexReader.close();
    }
    else std::cerr << "Unable to open FileIndex file";

    string XMLFileStr = "./Envi" + std::to_string(FileIndex) + ".xml";
    const char* XMLFile = XMLFileStr.c_str();    // Here we must give abstract path to the file
    if(!Backend.LoadAndInitSim(XMLFile))
    {
      std::cerr<<"World XML file path does not exist!"<<endl;
    }
    const int NumberOfTerrains = world.terrains.size();
    for (int i = 6; i < SimRobotObj.q.size(); i++)
    {
      int Linkid = world.RobotLinkID(0, i);
      Meshing::TriMesh LinkiTriMesh = world.GetGeometry(Linkid)->AsTriangleMesh();
      CollisionMesh LinkiTriMeshTopology(LinkiTriMesh);
      LinkiTriMeshTopology.InitCollisions();
      LinkiTriMeshTopology.CalcTriNeighbors();
      for (int j = 0; j < NumberOfTerrains; j++)
      {
        Meshing::TriMesh EnviTriMesh  = world.terrains[j]->geometry->AsTriangleMesh();
        CollisionMesh EnviTriMeshTopology(EnviTriMesh);
        EnviTriMeshTopology.InitCollisions();
        EnviTriMeshTopology.CalcTriNeighbors();
        CollisionMeshQuery CollisionMeshes(LinkiTriMeshTopology, EnviTriMeshTopology);
        bool CollsiionCheck = CollisionMeshes.Collide();
        if (CollsiionCheck == true)
        {
          std::printf("Link %d and Terrain %d are in collision!\n", i, j);
        }
      }
    }
    int a = 1;
  }
  return InitConfig;
}
