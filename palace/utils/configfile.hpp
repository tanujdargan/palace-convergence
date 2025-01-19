// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_CONFIG_FILE_HPP
#define PALACE_UTILS_CONFIG_FILE_HPP

#include <array>
#include <map>
#include <string>
#include <vector>
#include <nlohmann/json_fwd.hpp>

namespace palace::config
{

using json = nlohmann::json;

//
// Data structures for storing configuration file data.
//

namespace internal
{

template <typename DataType>
struct DataVector
{
protected:
  std::vector<DataType> vecdata = {};

public:
  [[nodiscard]] const auto &operator[](int i) const { return vecdata[i]; }
  [[nodiscard]] auto &operator[](int i) { return vecdata[i]; }
  [[nodiscard]] const auto &at(int i) const { return vecdata.at(i); }
  [[nodiscard]] auto &at(int i) { return vecdata.at(i); }
  [[nodiscard]] auto size() const { return vecdata.size(); }
  [[nodiscard]] auto empty() const { return vecdata.empty(); }
  [[nodiscard]] auto begin() const { return vecdata.begin(); }
  [[nodiscard]] auto end() const { return vecdata.end(); }
  [[nodiscard]] auto begin() { return vecdata.begin(); }
  [[nodiscard]] auto end() { return vecdata.end(); }
  [[nodiscard]] auto front() const { return vecdata.front(); }
  [[nodiscard]] auto back() const { return vecdata.back(); }
  [[nodiscard]] auto front() { return vecdata.front(); }
  [[nodiscard]] auto back() { return vecdata.back(); }
};

template <typename DataType>
struct DataMap
{
protected:
  // Map keys are the object indices for postprocessing.
  std::map<int, DataType> mapdata = {};

public:
  [[nodiscard]] const auto &operator[](int i) const { return mapdata[i]; }
  [[nodiscard]] auto &operator[](int i) { return mapdata[i]; }
  [[nodiscard]] const auto &at(int i) const { return mapdata.at(i); }
  [[nodiscard]] auto &at(int i) { return mapdata.at(i); }
  [[nodiscard]] auto size() const { return mapdata.size(); }
  [[nodiscard]] auto empty() const { return mapdata.empty(); }
  [[nodiscard]] auto begin() const { return mapdata.begin(); }
  [[nodiscard]] auto end() const { return mapdata.end(); }
  [[nodiscard]] auto begin() { return mapdata.begin(); }
  [[nodiscard]] auto end() { return mapdata.end(); }
};

// An ElementData consists of a list of attributes making up a single element of a
// potentially multielement boundary, and a direction and/or a normal defining the incident
// field. These are used for lumped ports, terminals, surface currents, and other boundary
// postprocessing objects.
struct ElementData
{
  // Vector defining the direction for this port. In a Cartesian system, "X", "Y", and "Z"
  // map to (1,0,0), (0,1,0), and (0,0,1), respectively.
  std::array<double, 3> direction{{0.0, 0.0, 0.0}};

  // Coordinate system that the normal vector is expressed in.
  enum class CoordinateSystem
  {
    CARTESIAN,
    CYLINDRICAL
  };
  CoordinateSystem coordinate_system = CoordinateSystem::CARTESIAN;

  // List of boundary attributes for this element.
  std::vector<int> attributes = {};
};

}  // namespace internal

struct ProblemData
{
public:
  // Simulation type.
  enum class Type
  {
    DRIVEN,
    EIGENMODE,
    ELECTROSTATIC,
    MAGNETOSTATIC,
    TRANSIENT
  };
  Type type = Type::DRIVEN;

  // Level of printing.
  int verbose = 1;

  // Output path for storing results.
  std::string output = "";

  void SetUp(json &config);
};

struct BoxRefinementData
{
  // Refinement levels.
  int ref_levels = 0;

  // Region bounding box limits [m].
  std::array<double, 3> bbmin{{0.0, 0.0, 0.0}}, bbmax{{0.0, 0.0, 0.0}};
};

struct SphereRefinementData
{
  // Refinement levels.
  int ref_levels = 0;

  // Sphere radius [m].
  double r = 0.0;

  // Sphere center [m].
  std::array<double, 3> center{{0.0, 0.0, 0.0}};
};

struct RefinementData
{
public:
  // Non-dimensional tolerance used to specify convergence of adaptive mesh refinement.
  double tol = 1.0e-2;

  // Maximum number of iterations to perform during adaptive mesh refinement.
  int max_it = 0;

  // If a refinement results in a greater number of DOFs than this value, no future
  // refinement will be allowed.
  int max_size = 0;

  // Whether or not to perform nonconformal adaptation.
  bool nonconformal = true;

  // Maximum difference in nonconformal refinements between two adjacent elements. Zero
  // implies there is no constraint on local nonconformity.
  int max_nc_levels = 1;

  // Dörfler update fraction. The set of marked elements is the minimum set that contains
  // update_fraction of the total error.
  double update_fraction = 0.7;

  // Maximum allowable ratio of number of elements across processors before rebalancing is
  // performed.
  double maximum_imbalance = 1.1;

  // Whether to save off results of each adaptation iteration as a subfolder within the post
  // processing directory.
  bool save_adapt_iterations = true;

  // Whether to write a (serial) mesh to file after mesh modification during AMR.
  bool save_adapt_mesh = false;

  // Parallel uniform mesh refinement levels.
  int uniform_ref_levels = 0;

  // Serial uniform mesh refinement levels.
  int ser_uniform_ref_levels = 0;

private:
  // Refinement data for mesh regions.
  std::vector<BoxRefinementData> box_list = {};
  std::vector<SphereRefinementData> sphere_list = {};

public:
  auto &GetBox(int i) { return box_list[i]; }
  const auto &GetBoxes() const { return box_list; }
  auto &GetBoxes() { return box_list; }

  auto &GetSphere(int i) { return sphere_list[i]; }
  const auto &GetSpheres() const { return sphere_list; }
  auto &GetSpheres() { return sphere_list; }

  void SetUp(json &model);
};

struct ModelData
{
public:
  // Mesh file.
  std::string mesh = "";

  // Mesh length unit and optional characteristic length scale for nondimensionalization
  // [m].
  double L0 = 1.0e-6;
  double Lc = -1.0;

  // Remove high-order curvature information from the mesh.
  bool remove_curvature = false;

  // Convert mesh to simplex elements.
  bool make_simplex = false;

  // Convert mesh to hexahedral elements (using tet-to-hex algorithm).
  bool make_hex = false;

  // Reorder elements based on spatial location after loading the serial mesh, which can
  // potentially increase memory coherency.
  bool reorder_elements = false;

  // Remove elements (along with any associated unattached boundary elements) from the mesh
  // which do not have any material properties specified.
  bool clean_unused_elements = true;

  // Split, or "crack", boundary elements lying on internal boundaries to decouple the
  // elements on either side.
  bool crack_bdr_elements = true;

  // When required, refine elements neighboring a split or crack in order to enable the
  // decoupling.
  bool refine_crack_elements = true;

  // Factor for displacing duplicated interior boundary elements, usually added just for
  // visualization.
  double crack_displ_factor = 1.0e-3;

  // Add new boundary elements for faces are on the computational domain boundary or which
  // have attached elements on either side with different domain attributes.
  bool add_bdr_elements = true;

  // Call MFEM's ReorientTetMesh as a check of mesh orientation after partitioning.
  bool reorient_tet_mesh = false;

  // Partitioning file (if specified, does not compute a new partitioning).
  std::string partitioning = "";

  // Object controlling mesh refinement.
  RefinementData refinement = {};

  void SetUp(json &config);
};

// Store symmetric matrix data as set of outer products: Σᵢ sᵢ * vᵢ *  vᵢᵀ.
template <std::size_t N>
struct SymmetricMatrixData
{
public:
  std::array<double, N> s;
  std::array<std::array<double, N>, N> v;

  SymmetricMatrixData(double diag)
  {
    s.fill(diag);
    std::size_t i = 0;
    for (auto &x : v)
    {
      x.fill(0.0);
      x[i++] = 1.0;
    }
  }
};

struct MaterialData
{
public:
  // Relative permeability.
  SymmetricMatrixData<3> mu_r = 1.0;

  // Relative permittivity.
  SymmetricMatrixData<3> epsilon_r = 1.0;

  // Loss tangent.
  SymmetricMatrixData<3> tandelta = 0.0;

  // Conductivity [S/m].
  SymmetricMatrixData<3> sigma = 0.0;

  // London penetration depth [m].
  double lambda_L = 0.0;

  // List of domain attributes for this material.
  std::vector<int> attributes = {};
};

struct DomainMaterialData : public internal::DataVector<MaterialData>
{
public:
  void SetUp(json &domains);
};

struct DomainEnergyData
{
public:
  // List of domain attributes for this domain postprocessing index.
  std::vector<int> attributes = {};
};

struct DomainEnergyPostData : public internal::DataMap<DomainEnergyData>
{
public:
  void SetUp(json &postpro);
};

struct ProbeData
{
public:
  // Physical space coordinates for the probe location [m].
  std::array<double, 3> center{{0.0, 0.0, 0.0}};
};

struct ProbePostData : public internal::DataMap<ProbeData>
{
public:
  void SetUp(json &postpro);
};

struct DomainPostData
{
public:
  // List of all postprocessing domain attributes.
  std::vector<int> attributes = {};

  // Domain postprocessing objects.
  DomainEnergyPostData energy;
  ProbePostData probe;

  void SetUp(json &domains);
};

struct DomainData
{
public:
  // List of all domain attributes (excluding postprocessing).
  std::vector<int> attributes = {};

  // Domain objects.
  DomainMaterialData materials = {};
  DomainPostData postpro = {};

  void SetUp(json &config);
};

struct PecBoundaryData
{
public:
  // List of boundary attributes with PEC boundary conditions.
  std::vector<int> attributes = {};

  [[nodiscard]] auto empty() const { return attributes.empty(); }

  void SetUp(json &boundaries);
};

struct PmcBoundaryData
{
public:
  // List of boundary attributes with PMC boundary conditions.
  std::vector<int> attributes = {};

  [[nodiscard]] auto empty() const { return attributes.empty(); }

  void SetUp(json &boundaries);
};

struct WavePortPecBoundaryData
{
public:
  // List of boundary attributes with PEC boundary conditions for wave ports.
  std::vector<int> attributes = {};

  [[nodiscard]] auto empty() const { return attributes.empty(); }

  void SetUp(json &boundaries);
};

struct FarfieldBoundaryData
{
public:
  // Approximation order for farfield ABC.
  int order = 1;

  // List of boundary attributes with farfield absortbing boundary conditions.
  std::vector<int> attributes = {};

  [[nodiscard]] auto empty() const { return attributes.empty(); }

  void SetUp(json &boundaries);
};

struct ConductivityData
{
public:
  // Electrical conductivity of the conductor [S/m].
  double sigma = 0.0;

  // Conductor relative permeability.
  double mu_r = 1.0;

  // Optional conductor thickness [m].
  double h = 0.0;

  // Optional flag for an external boundary surface, relevant for the thickness correction.
  bool external = false;

  // List of boundary attributes for this surface conductivity boundary condition.
  std::vector<int> attributes = {};
};

struct ConductivityBoundaryData : public internal::DataVector<ConductivityData>
{
public:
  void SetUp(json &boundaries);
};

struct ImpedanceData
{
public:
  // Boundary surface resistance, inductance, and capacitance [Ω/sq, H/sq, F/sq].
  double Rs = 0.0;
  double Ls = 0.0;
  double Cs = 0.0;

  // List of boundary attributes for this impedance boundary condition.
  std::vector<int> attributes = {};
};

struct ImpedanceBoundaryData : public internal::DataVector<ImpedanceData>
{
public:
  void SetUp(json &boundaries);
};

struct LumpedPortData
{
public:
  // Port circuit resistance, inductance, and capacitance [Ω/sq, H/sq, F/sq].
  double R = 0.0;
  double L = 0.0;
  double C = 0.0;

  // Port surface resistance, inductance, and capacitance [Ω/sq, H/sq, F/sq].
  double Rs = 0.0;
  double Ls = 0.0;
  double Cs = 0.0;

  // Flag for source term in driven and transient simulations.
  bool excitation = false;

  // Flag for boundary damping term in driven and transient simulations.
  bool active = true;

  // For each lumped port index, each element contains a list of attributes making up a
  // single element of a potentially multielement lumped port.
  std::vector<internal::ElementData> elements = {};
};

struct LumpedPortBoundaryData : public internal::DataMap<LumpedPortData>
{
public:
  void SetUp(json &boundaries);
};

struct PeriodicData
{
public:
  // Vector defining the direction and distance for this periodic boundary condition.
  std::array<double, 3> translation = {0.0, 0.0, 0.0};
  // List of boundary donor attributes for this periodic boundary condition.
  std::vector<int> donor_attributes = {};
  // List of boundary receiver attributes for this periodic boundary condition.
  std::vector<int> receiver_attributes = {};
};

struct PeriodicBoundaryData : public internal::DataVector<PeriodicData>
{
public:
  void SetUp(json &boundaries);
};

struct WavePortData
{
public:
  // Mode index for the numeric wave port.
  int mode_idx = 1;

  // Port offset for de-embedding [m].
  double d_offset = 0.0;

  // Eigenvalue solver type for boundary mode calculation.
  enum class EigenSolverType
  {
    DEFAULT,
    SLEPC,
    ARPACK
  };
  EigenSolverType eigen_type = EigenSolverType::DEFAULT;

  // Flag for source term in driven and transient simulations.
  bool excitation = false;

  // Flag for boundary damping term in driven and transient simulations.
  bool active = true;

  // List of boundary attributes for this wave port.
  std::vector<int> attributes = {};

  // Maximum number of iterations in linear solver.
  int ksp_max_its = 45;

  // Tolerance for linear solver.
  double ksp_tol = 1e-8;

  // Tolerance for eigenvalue solver.
  double eig_tol = 1e-6;

  // Print level for linear and eigenvalue solvers.
  int verbose = 0;
};

struct WavePortBoundaryData : public internal::DataMap<WavePortData>
{
public:
  void SetUp(json &boundaries);
};

struct SurfaceCurrentData
{
public:
  // For each surface current source index, each element contains a list of attributes
  // making up a single element of a potentially multielement current source.
  std::vector<internal::ElementData> elements = {};
};

struct SurfaceCurrentBoundaryData : public internal::DataMap<SurfaceCurrentData>
{
public:
  void SetUp(json &boundaries);
};

struct SurfaceFluxData
{
public:
  // Surface flux type.
  enum class Type
  {
    ELECTRIC,
    MAGNETIC,
    POWER
  };
  Type type = Type::ELECTRIC;

  // Flag for whether or not to consider the boundary as an infinitely thin two-sided
  // boundary for postprocessing.
  bool two_sided = false;

  // Coordinates of a point away from which to compute the outward flux (for orienting the
  // surface normal) [m].
  std::array<double, 3> center{{0.0, 0.0, 0.0}};

  // Flag which indicates whether or not the center point was specified.
  bool no_center = true;

  // List of boundary attributes for this surface flux postprocessing index.
  std::vector<int> attributes = {};
};

struct SurfaceFluxPostData : public internal::DataMap<SurfaceFluxData>
{
public:
  void SetUp(json &postpro);
};

struct InterfaceDielectricData
{
public:
  // Type of interface dielectric for computing electric field energy participation ratios.
  enum class Type
  {
    DEFAULT,
    MA,
    MS,
    SA
  };
  Type type = Type::DEFAULT;

  // Dielectric interface thickness [m].
  double t = 0.0;

  // Relative permittivity.
  double epsilon_r = 0.0;

  // Loss tangent.
  double tandelta = 0.0;

  // List of boundary attributes for this interface dielectric postprocessing index.
  std::vector<int> attributes = {};
};

struct InterfaceDielectricPostData : public internal::DataMap<InterfaceDielectricData>
{
public:
  void SetUp(json &postpro);
};

struct BoundaryPostData
{
public:
  // List of all postprocessing boundary attributes.
  std::vector<int> attributes = {};

  // Boundary postprocessing objects.
  SurfaceFluxPostData flux = {};
  InterfaceDielectricPostData dielectric = {};

  void SetUp(json &boundaries);
};

struct BoundaryData
{
public:
  // List of all boundary attributes (excluding postprocessing).
  std::vector<int> attributes = {};

  // Boundary objects.
  PecBoundaryData pec = {};
  PmcBoundaryData pmc = {};
  WavePortPecBoundaryData auxpec = {};
  FarfieldBoundaryData farfield = {};
  ConductivityBoundaryData conductivity = {};
  ImpedanceBoundaryData impedance = {};
  LumpedPortBoundaryData lumpedport = {};
  WavePortBoundaryData waveport = {};
  SurfaceCurrentBoundaryData current = {};
  PeriodicBoundaryData periodic = {};
  BoundaryPostData postpro = {};

  void SetUp(json &config);
};

struct DrivenSolverData
{
public:
  // Lower bound of frequency sweep [GHz].
  double min_f = 0.0;

  // Upper bound of frequency sweep [GHz].
  double max_f = 0.0;

  // Step size for frequency sweep [GHz].
  double delta_f = 0.0;

  // Step increment for saving fields to disk.
  int delta_post = 0;

  // Restart iteration for a partial sweep.
  int rst = 1;

  // Error tolerance for enabling adaptive frequency sweep.
  double adaptive_tol = 0.0;

  // Maximum number of frequency samples for adaptive frequency sweep.
  int adaptive_max_size = 0;

  // Memory required for adaptive sampling convergence.
  int adaptive_memory = 2;

  void SetUp(json &solver);
};

struct JosephsonElementData
{
    // Location of the Josephson element in 3D space
    std::vector<double> location;
    // Unique identifier for the element
    int index;
    // Material properties
    double critical_current;  // Critical current in amperes
    double junction_capacitance;  // Junction capacitance in farads
};

struct EigenmodeSolverData
{
public:
    // Number of eigenvalues to compute
    int n = 1;
    // Number of fields to write to disk
    int n_post = 0;
    // Convergence tolerance
    double tol = 1.0e-6;
    // Target frequency [GHz]
    double target = 0.0;
    // Use a constant vector as initial guess
    bool start_vector_constant = false;
    // Josephson element configuration
    std::vector<JosephsonElementData> josephson_elements;
    // Field convergence threshold for Josephson elements
    double field_convergence_threshold = 1e-6;
    // Maximum number of iterations for field convergence
    int max_field_iterations = 100;
    
    void SetUp(json &solver);
};

struct JosephsonConfig 
{
    // Vector of locations to probe field values
    std::vector<mfem::Vector> locations;
    // Convergence threshold for field values
    double convergence_threshold = 1e-6;
    // Maximum iterations for field convergence
    int max_iterations = 100;
};

struct SolverData
{
public:
  // Approximation order.
  int order = 1;

  // Order above which to use partial assembly instead of full assembly.
  int pa_order_threshold = 1;

  // Include the order of det(J) in the order of accuracy for quadrature rule selection.
  bool q_order_jac = false;

  // Additional quadrature order of accuracy (in addition to 2p or 2p + order(|J|)) for
  // quadrature rule selection.
  int q_order_extra = 0;

  // Device used to configure MFEM.
  enum class Device
  {
    CPU,
    GPU,
    DEBUG
  };
  Device device = Device::CPU;

  // Backend for libCEED (https://libceed.org/en/latest/gettingstarted/#backends).
  std::string ceed_backend = "";

  // Solver objects.
  DrivenSolverData driven = {};
  EigenmodeSolverData eigenmode = {};
  ElectrostaticSolverData electrostatic = {};
  MagnetostaticSolverData magnetostatic = {};
  TransientSolverData transient = {};
  LinearSolverData linear = {};

  void SetUp(json &config);
};

}  // namespace palace::config

#endif  // PALACE_UTILS_CONFIGFILE_HPP
