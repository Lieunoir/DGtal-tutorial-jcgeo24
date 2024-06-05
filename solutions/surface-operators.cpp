#include <DGtal/base/Common.h>
#include <DGtal/dec/InterpolatedCorrectedCalculus.h>
#include <DGtal/dec/NormalCorrectedFEM.h>
#include <DGtal/helpers/Shortcuts.h>
#include <DGtal/helpers/ShortcutsGeometry.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/math/linalg/DirichletConditions.h>
#include <DGtal/math/linalg/EigenSupport.h>
#include <iostream>

#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

using namespace DGtal;
using namespace Z3i;

// Using standard 3D digital space.
typedef Shortcuts<KSpace> SH3;
typedef ShortcutsGeometry<KSpace> SHG3;
typedef SurfaceMesh<RealPoint, RealVector> SurfMesh;
typedef SurfMesh::Face Face;
typedef SurfMesh::Vertex Vertex;
typedef SurfMesh::Index Index;
typedef NormalCorrectedFEM<EigenLinearAlgebraBackend, SH3::RealPoint,
                           SH3::RealVector>
    ncFEM;
typedef InterpolatedCorrectedCalculus<EigenLinearAlgebraBackend, SH3::RealPoint,
                                      SH3::RealVector>
    CC;
typedef DirichletConditions<EigenLinearAlgebraBackend> DC;
typedef EigenLinearAlgebraBackend::DenseVector DenseVector;

// Global variables to make easier GUI stuff.
polyscope::SurfaceMesh *psMesh;
SurfMesh surfmesh;
float GridStep = 0.5;

/// Create an implicit shape \a polynomial digitized at gridstep \a h
/// @param polynomial the implicit function as a  multivariate polynomial
/// string.
/// @param h the chosen digitization gridstep
void createShape(std::string polynomial, double h) {
  auto params = SH3::defaultParameters() | SHG3::defaultParameters() |
                SHG3::parametersGeometryEstimation();
  params("surfaceComponents", "All");
  params("polynomial", polynomial)("minAABB", -10.0)("maxAABB", 10.0)(
      "offset", 1.0)("gridstep", h);
  auto shape = SH3::makeImplicitShape3D(params);
  auto dshape = SH3::makeDigitizedImplicitShape3D(shape, params);
  auto K = SH3::getKSpace(params);
  auto binary_image = SH3::makeBinaryImage(dshape, params);
  auto surface = SH3::makeDigitalSurface(binary_image, K, params);
  auto primalSurface = SH3::makePrimalSurfaceMesh(surface);
  auto surfels = SH3::getSurfelRange(surface, params);
  auto true_normals = SHG3::getNormalVectors(shape, K, surfels, params);

  // Need to convert the faces
  std::vector<std::vector<SH3::SurfaceMesh::Vertex>> faces;
  std::vector<RealPoint> positions;
  std::vector<RealPoint> smooth_positions;

  for (auto face = 0; face < primalSurface->nbFaces(); ++face)
    faces.push_back(primalSurface->incidentVertices(face));

  // Embed lattice points according to gridstep.
  positions = primalSurface->positions();
  for (auto &x : positions)
    x *= h;

  // Create DGtal surface mesh object.
  surfmesh =
      SurfMesh(positions.begin(), positions.end(), faces.begin(), faces.end());
  surfmesh.setFaceNormals(true_normals.begin(), true_normals.end());
  std::cout << surfmesh << std::endl;
  std::cout << "number of non-manifold Edges = "
            << surfmesh.computeNonManifoldEdges().size() << std::endl;
  // Create rendered polyscope surface.
  psMesh = polyscope::registerSurfaceMesh("digital surface", positions, faces);
  psMesh->addFaceVectorQuantity("True normal vector field", true_normals);
}

void estimateAreaMassmatrix() {
  ncFEM calculus(surfmesh);
  ncFEM::LinearOperator M = calculus.M0();
  DenseVector ones = DenseVector::Ones(surfmesh.nbVertices());
  std::cout << "Estimated area : " << ones.transpose() * M * ones << std::endl;
}

void diffuseHeat() {
  ncFEM calculus(surfmesh);
  ncFEM::LinearOperator M = calculus.M0();
  ncFEM::LinearOperator L = calculus.L0();
  DenseVector source = DenseVector::Zero(surfmesh.nbVertices());
  source(1) = 1.;
  EigenLinearAlgebraBackend::SolverSimplicialLDLT solver;
  double dt = GridStep * GridStep;
  solver.compute(M - dt * L);
  DenseVector diffused = solver.solve(M * source);
  psMesh->addVertexScalarQuantity("Source", source);
  psMesh->addVertexScalarQuantity("Diffused heat", diffused);
}

void smoothInterpolation() {
  ncFEM calculus(surfmesh);
  ncFEM::LinearOperator L = calculus.L0();
  DenseVector g = DenseVector::Zero(surfmesh.nbVertices());
  DC::IntegerVector b = DC::IntegerVector::Zero(g.rows());

  // Add 4 arbitrary sources
  g(1) = 10.;
  g(50) = -10.;
  g(100) = 10.;
  g(150) = -10.;
  b(1) = 1;
  b(50) = 1;
  b(100) = 1;
  b(150) = 1;

  // Solve Δu=0 with g as boundary conditions
  ncFEM::LinearAlgebraBackend::SolverSimplicialLDLT solver;
  ncFEM::LinearOperator L_dirichlet = DC::dirichletOperator(L, b);
  solver.compute(L_dirichlet);
  ASSERT(solver.info() == Eigen::Success);
  DenseVector g_dirichlet = DC::dirichletVector(L, g, b, g);
  DenseVector x_dirichlet = solver.solve(g_dirichlet);
  DenseVector u = DC::dirichletSolution(x_dirichlet, b, g);
  psMesh->addVertexScalarQuantity("To interpolate", g);
  psMesh->addVertexScalarQuantity("Interpolated", u);
}

void heatGeodesics() {
  surfmesh.computeVertexNormalsFromFaceNormals();
  CC calculus(surfmesh);
  CC::LinearOperator L = calculus.L0();
  CC::LinearOperator M = calculus.M0();
  CC::LinearOperator Grad = calculus.Sharp() * calculus.D0();
  CC::LinearOperator Div =
      calculus.D0().transpose() * calculus.M1() * calculus.Flat();

  DenseVector source = DenseVector::Zero(surfmesh.nbVertices());
  source(0) = 1.;
  EigenLinearAlgebraBackend::SolverSimplicialLDLT solver;
  double dt = GridStep * GridStep;
  solver.compute(M - dt * L);
  EigenLinearAlgebraBackend::DenseVector diffused = solver.solve(M * source);
  psMesh->addVertexScalarQuantity("Geodesics heat source", source);
  psMesh->addVertexScalarQuantity("Geodesics diffused heat", diffused);

  DenseVector grads = Grad * diffused;
  for (int i = 0; i < surfmesh.nbFaces(); i++) {
    double c0 = grads(3 * i);
    double c1 = grads(3 * i + 1);
    double c2 = grads(3 * i + 2);
    double norm = sqrt(c0 * c0 + c1 * c1 + c2 * c2);
    grads(3 * i) = c0 / norm;
    grads(3 * i + 1) = c1 / norm;
    grads(3 * i + 2) = c2 / norm;
  }
  DenseVector divs = Div * grads;

  DenseVector g = DenseVector::Zero(surfmesh.nbVertices());
  DC::IntegerVector b = DC::IntegerVector::Zero(g.rows());
  g(0) = 0.;
  b(0) = 1;
  CC::LinearOperator L_dirichlet = DC::dirichletOperator(L, b);
  solver.compute(L_dirichlet);
  ASSERT(solver.info() == Eigen::Success);
  DenseVector g_dirichlet = DC::dirichletVector(L, divs, b, g);
  DenseVector x_dirichlet = solver.solve(g_dirichlet);
  DenseVector u = DC::dirichletSolution(x_dirichlet, b, g);
  psMesh->addVertexScalarQuantity("Geodesics", u);
}

/// Defines the GUI buttons and reactions.
void myCallback() {
  if (ImGui::Button("Sphere"))
    createShape("sphere9", GridStep);
  ImGui::SameLine();
  if (ImGui::Button("Torus"))
    createShape("torus", GridStep);
  ImGui::SameLine();
  if (ImGui::Button("Cylinder"))
    createShape("x^2-2*x*y+y^2+z^2-25", GridStep);
  ImGui::SliderFloat("Gridstep h parameter", &GridStep, 0.025, 2.0);
  if (ImGui::Button("Estimate Area"))
    estimateAreaMassmatrix();
  if (ImGui::Button("Heat diffusion"))
    diffuseHeat();
  if (ImGui::Button("Smooth interpolation"))
    smoothInterpolation();
  if (ImGui::Button("Heat geodesics"))
    heatGeodesics();
}

int main() {
  // Gives you the list of predefined shapes
  auto L = SH3::getPolynomialList();
  for (const auto &e : L)
    std::cout << e.first << " : " << e.second << std::endl;

  // Initialize polyscope
  polyscope::init();

  // Create shape
  createShape("sphere9", GridStep);

  // Set the callback function
  polyscope::state::userCallback = myCallback;
  polyscope::show();
  return EXIT_SUCCESS;
}
