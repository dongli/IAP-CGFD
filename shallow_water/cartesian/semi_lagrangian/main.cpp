/**
 *  Description:
 *
 *    This is an implementation of a two-time-level semi-implicit semi-Lagrangian
 *    shallow-water model in Cartesian geometry described in Staniforth and Cote (1991).
 *
 *  References:
 *
 *    - Andrew Staniforth, and Jean Cote, 1991: Semi-Lagrangian Integration Schemes for
 *      Atmospheric Models - A Review. Monthly Weather Review, 119, 2206-2223.
 */

#include "geomtk/Cartesian.h"

using namespace std;

int main(int argc, const char *argv[])
{
    ConfigManager configManager;
    Domain domain(2);
    Mesh mesh(domain, 2);
    Field<double, 3> u, v, h;
    Field<arma::vec::fixed<2> > a;
    TimeManager timeManager;
    IOManager io;
    int outputFileIdx;
    TimeLevelIndex<3> lp1, lp05, l0, lm1;

    Regrid regrid(mesh);
    RegridMethod velocityRegridMethod;

    double dt, dx, dy;
    double g;
    int numIter;
    string outputPattern = "semi_lagrangian.%3s.nc";
    RegridMethod velocityInterpMethod;
    RegridMethod advectInterpMethod;

    if (argc != 2) {
        REPORT_ERROR("Configure file is needed!");
    }
    configManager.parse(argv[1]);
    dt = configManager.getValue("semi_lagrangian", "dt", 1);
    dx = configManager.getValue("semi_lagrangian", "dx", 0.01);
    dy = configManager.getValue("semi_lagrangian", "dy", 0.01);
    g = configManager.getValue("semi_lagrangian", "gravitational_accelaration", 9.8);
    numIter = configManager.getValue("semi_lagrangian", "num_iteration", 2);
    outputPattern = configManager.getValue("semi_lagrangian", "output_pattern", outputPattern);
    velocityInterpMethod = Regrid::methodFromString(configManager.getValue<string>("semi_lagrangian", "velocity_interp_method"));
    advectInterpMethod = Regrid::methodFromString(configManager.getValue<string>("semi_lagrangian", "advect_interp_method"));

    // Set the one dimensional space axis.
    domain.setAxis(0, "x", "x axis", "m",
                   0, geomtk::BndType::PERIODIC,
                   1, geomtk::BndType::PERIODIC);
    domain.setAxis(1, "y", "y axis", "m",
                   0, geomtk::BndType::PERIODIC,
                   1, geomtk::BndType::PERIODIC);

    // Set the discrete mesh on the domain.
    mesh.init(domain.axisSpan(0)/dx, domain.axisSpan(1)/dy);

    // Set the time manager.
    Time startTime(0*geomtk::TimeUnit::SECONDS);
    Time endTime(200*geomtk::TimeUnit::SECONDS);
    timeManager.init(startTime, endTime, dt);

    // Set up velocity and density fields.
    u.create("u", "m s-1", "velocity component along x axis", mesh, CENTER, 2, true);
    v.create("v", "m s-1", "velocity component along y axis", mesh, CENTER, 2, true);
    h.create("h", "m", "height", mesh, CENTER, 2);
    a.create("alpha", "m", "displacement", mesh, CENTER, 2);

    // Set up initial guess for displacements.
    for (int j = mesh.js(FULL); j < mesh.je(FULL); ++j) {
        for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
            a(i, j).zeros();
        }
    }

    // Set up initial conditions.
    l0 = lm1+1;
    for (int j = mesh.js(FULL); j < mesh.je(FULL); ++j) {
        for (int i = mesh.is(FULL); i < mesh.ie(FULL); ++i) {
            u(l0,  i, j) = 0; v(l0,  i, j) = 0;
            u(lm1, i, j) = 0; v(lm1, i, j) = 0;
            
        }
    }

    // Set up IO manager.
    io.init(timeManager);
    outputFileIdx = io.registerOutputFile(mesh, outputPattern, geomtk::TimeStepUnit::STEP, 1);
    io.registerField(outputFileIdx, "double", FULL_DIMENSION, {&u, &v, &h});
    io.output<double, 3>(outputFileIdx, l0, {&u, &v, &h});

    // Run the main loop.
    while (!timeManager.isFinished()) {
        lp1 = l0+1; lp05 = l0+0.5; lm1 = l0-1;
        // Extrapolate velocity on t+1/2*dt time level.
        for (int j = mesh.js(FULL); j < mesh.je(FULL); ++j) {
            for (int i = mesh.is(FULL); i < mesh.ie(FULL); ++i) {
                u(lp05, i, j) = 1.5*u(l0, i, j)-0.5*u(lm1, i, j);
                v(lp05, i, j) = 1.5*v(l0, i, j)-0.5*v(lm1, i, j);
            }
        }
        u.applyBndCond(lp05);
        v.applyBndCond(lp05);
        // Calculate displacements.
        for (int iter = 0; iter < numIter; ++iter) {
            for (int j = mesh.js(FULL); j < mesh.je(FULL); ++j) {
                for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
                    const SpaceCoord &x0 = mesh.gridCoord(CENTER, i, j);
                    SpaceCoord x(2);
                    x() = x0()-0.5*a(i, j);
                    MeshIndex meshIdx(2);
                    meshIdx.locate(mesh, x);
                    double U, V;
                    regrid.run(velocityRegridMethod, lp05, u, x, U, &meshIdx);
                    regrid.run(velocityRegridMethod, lp05, v, x, V, &meshIdx);
                    a(i, j)(0) = dt*U;
                    a(i, j)(1) = dt*V;
                }
            }
        }
        timeManager.advance(); l0.shift();
        io.output<double, 3>(outputFileIdx, l0, {&u, &v, &h});
    }

    return 0;
}
