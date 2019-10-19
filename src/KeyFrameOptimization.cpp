#include "NonlinearOptimizerInfo.h"
#include "CommonHeader.h"

// This function is used to generate key frames for walking gait.

static Robot SimRobotObj;
static std::vector<LinkInfo> RobotLinkInfo;
static std::vector<ContactStatusInfo> RobotContactInfo;
static SignedDistanceFieldInfo SDFInfo;
static std::vector<double> RobotConfigRef;
static std::vector<Vector3> ContactPosRef;

struct KeyFrameOpt: public NonlinearOptimizerInfo
{
  // This is the class struct for the optimizaiton of configuraiton variables

  KeyFrameOpt():NonlinearOptimizerInfo(){};

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
      std::vector<double> F_val = KeyFrameObjNCons(*n, *neF, x_vec);
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
  static std::vector<double> KeyFrameObjNCons(const int&nVar, const int &nObjNCons, const std::vector<double>& ConfigOpt)
  {
    // This funciton provides the constraint for the configuration variable

    std::vector<double> F(nObjNCons);
    Config ConfigOptNew(ConfigOpt);
    SimRobotObj.UpdateConfig(ConfigOptNew);     // Here both the SimRobot.q and robot frames have already been updated.

    std::vector<Vector3> ContactPosNow, ActVelocitiesRef;
    std::vector<Matrix> ActJacobiansRef;
    std::vector<double> ActDistsRef;
    std::vector<int> ActStatus = ActContactNJacobian(SimRobotObj, RobotLinkInfo, RobotContactInfo, ContactPosNow, ActVelocitiesRef, ActDistsRef, ActJacobiansRef, SDFInfo);

    // The objective is to make sure the contact points are not moving
    double ConfigVia = 0.0;
    for (int i = 0; i <ContactPosNow.size(); i++)
    {
      double ContactVia_i_x = ContactPosRef[i].x - ContactPosNow[i].x;
      double ContactVia_i_y = ContactPosRef[i].y - ContactPosNow[i].y;
      double ContactVia_i_z = ContactPosRef[i].z - ContactPosNow[i].z;
      double ContactVia_i = ContactVia_i_x * ContactVia_i_x + ContactVia_i_y * ContactVia_i_y + ContactVia_i_z * ContactVia_i_z;
      ConfigVia +=ContactVia_i;
    }
    // Frame3:
    // Frame5:

    Vector3 Frame3(0.36662071646319372, 0.26861018572382089, 0.69419804801816232);
    Vector3 Frame5(0.49720115523986069, 0.20236116877678717, 0.74778429631628218);

    Vector3 COMPos(0.0, 0.0, 0.0), COMVel(0.0, 0.0, 0.0), COMAcc(0.0, 0.0, 0.0);
    CentroidalState(SimRobotObj, COMPos, COMVel);

    double COMPos_x = COMPos.x - 0.5 * (Frame3.x + Frame5.x);
    double COMPos_y = COMPos.y - 0.5 * (Frame3.y + Frame5.y);
    double COMPos_z = COMPos.z - 0.5 * (Frame3.z + Frame5.z);
    double COMPos_Diff = COMPos_x * COMPos_x + COMPos_y * COMPos_y + COMPos_z * COMPos_z;
    ConfigVia+=COMPos_Diff;
    F[0] =  ConfigVia;
    // Make sure that active end effectors have zero relative signed distance.
    int ConstraintIndex = 1;
    for (int i = 0; i < RobotLinkInfo.size(); i++)
    {
      for (int j = 0; j < RobotLinkInfo[i].LocalContacts.size(); j++)
      {
        Vector3 LinkiPjPos;
        SimRobotObj.GetWorldPosition(RobotLinkInfo[i].LocalContacts[j], RobotLinkInfo[i].LinkIndex, LinkiPjPos);
        F[ConstraintIndex] = SDFInfo.SignedDistance(LinkiPjPos);      ConstraintIndex = ConstraintIndex + 1;
      }
    }

    // The last constraint is that all the links have to remain nonnegative with respect to the environment.
    // The robot link COM position
    for (int i = 0; i < nVar-6; i++)
    {
      Vector3 Link_i_COM;
      SimRobotObj.links[i+6].GetWorldCOM(Link_i_COM);
      F[ConstraintIndex] = SDFInfo.SignedDistance(Link_i_COM);     // Here 0.001 is the epsilon addition
      ConstraintIndex = ConstraintIndex + 1;
    }

    // The last constraint is to make sure that Center of Mass remains to be within SP
    std::vector<Vector3> SPVertices;
    for (int i = 0; i < RobotLinkInfo.size(); i++)
    {
      // if(i == 1)
      // {
      //   continue;
      // }

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
    F[ConstraintIndex] = SPObj.ProjPoint2EdgeDist(COM_Pos) - 0.1;
    ConstraintIndex = ConstraintIndex + 1;

    return F;
  }
};

static void KeyFrameOptFn(std::vector<double> &RobotConfig)
{
  // This function is used to optimize the robot configuraiton variable such that certain contact is active.
  KeyFrameOpt KeyFrameOptProblem;

  // Static Variable Substitution
  int n = SimRobotObj.q.size();
  int neF = 1;                                                                  // Cost function on the norm difference between the reference configuration and the optimized configuration
  for (int i = 0; i < RobotContactInfo.size(); i++)
  {
    for (int j = 0; j < RobotContactInfo[i].LocalContactStatus.size(); j++)    // Each contact has a distance constraint.
    {
      neF = neF + 1;
    }
  }
  neF = neF + n - 6;
  neF = neF + 1;      // Add one more constraint on the CoM within SP
  KeyFrameOptProblem.InnerVariableInitialize(n, neF);

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
  KeyFrameOptProblem.VariableBoundsUpdate(xlow_vec, xupp_vec);

  /*
    Initialize the bounds of variables
  */
  std::vector<double> Flow_vec(neF), Fupp_vec(neF);
  Flow_vec[0] = 0;
  Fupp_vec[0] = 1e20;           // The default idea is that the objective should always be nonnegative
  int ConstraintIndex = 1;
  for (int i = 0; i < RobotContactInfo.size(); i++)
  {
    for (int j = 0; j < RobotContactInfo[i].LocalContactStatus.size(); j++)
    {
      // Distance Constraint: scalar
      switch (RobotContactInfo[i].LocalContactStatus[j])
      {
        case 1:     // This means that the certain constraint is active
        Flow_vec[ConstraintIndex] = 0;
        Fupp_vec[ConstraintIndex] = 0;
        ConstraintIndex = ConstraintIndex + 1;
        break;
        case 0:
        Flow_vec[ConstraintIndex] = 0;
        Fupp_vec[ConstraintIndex] = 1e20;         // This is due to the nonpenetraition consideration.
        ConstraintIndex = ConstraintIndex + 1;
        break;
      }
    }
  }
  for (int i = 0; i < SimRobotObj.q.size()-6; i++)
  {
    Flow_vec[ConstraintIndex] = 0;
    Fupp_vec[ConstraintIndex] = 1e20;         // This is due to the nonpenetraition consideration.
    ConstraintIndex = ConstraintIndex + 1;
  }

  // Projected Center of Mass within Support Polygon
  Flow_vec[ConstraintIndex] = 0;
  Fupp_vec[ConstraintIndex] = 1e20;
  ConstraintIndex = ConstraintIndex + 1;

  KeyFrameOptProblem.ConstraintBoundsUpdate(Flow_vec, Fupp_vec);

  /*
    Initialize the seed guess
  */
  KeyFrameOptProblem.SeedGuessUpdate(RobotConfigRef);

  /*
  Given a name of this problem for the output
  */
  KeyFrameOptProblem.ProblemNameUpdate("KeyFrameOptProblem", 0);

  KeyFrameOptProblem.NonlinearProb.setIntParameter("Iterations limit", 1000000);
  KeyFrameOptProblem.NonlinearProb.setIntParameter("Major iterations limit", 250);
  KeyFrameOptProblem.NonlinearProb.setIntParameter("Major print level", 0);
  KeyFrameOptProblem.NonlinearProb.setIntParameter("Minor print level", 0);
  /*
  ProblemOptions seting

  */
  // Solve with Finite-Difference
  KeyFrameOptProblem.ProblemOptionsUpdate(0, 3);
  KeyFrameOptProblem.Solve(RobotConfig);
  return;
}

bool KeyFrameOptimization(Robot& _SimRobotObj, const std::vector<LinkInfo> & _RobotLinkInfo, const std::vector<ContactStatusInfo> &  _RobotContactInfo, const SignedDistanceFieldInfo& _SDFInfo, const std::vector<double>& _RobotConfigRef, std::vector<double> & RobotConfig)
{
  SimRobotObj = _SimRobotObj;
  RobotLinkInfo = _RobotLinkInfo;
  RobotContactInfo = _RobotContactInfo;
  SDFInfo = _SDFInfo;
  RobotConfigRef = _RobotConfigRef;

  std::vector<Vector3> ActContactPositionsRef, ActVelocitiesRef;
  std::vector<Matrix> ActJacobiansRef;
  std::vector<double> ActDistsRef;
  std::vector<int> ActStatus = ActContactNJacobian(_SimRobotObj, RobotLinkInfo, RobotContactInfo, ActContactPositionsRef, ActVelocitiesRef, ActDistsRef, ActJacobiansRef, SDFInfo);

  ContactPosRef = ActContactPositionsRef;

  KeyFrameOptFn(RobotConfig);

  return true;
}

std::vector<double> KeyFrameMirror(const std::vector<double> &RobotConfig)
{
  // This function is used to take in the current robot's configuration and output a mirror robot configuration.
  std::vector<double> RobotConfigMir = RobotConfig;
  for (int i = 6; i < 11; i++)
  {
    RobotConfigMir[i] = RobotConfig[6 + i];
  }
  for (int i = 12; i < 17; i++)
  {
    RobotConfigMir[i] = RobotConfig[i-6];
  }
  for (int i = 22; i < 28; i++)
  {
    RobotConfigMir[i] = RobotConfig[i+7];
  }
  for (int i = 29; i < 35; i++)
  {
    RobotConfigMir[i] = RobotConfig[i-7];
  }
  return RobotConfigMir;
}
