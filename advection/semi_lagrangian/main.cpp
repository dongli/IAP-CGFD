#include "geomtk/Cartesian.h"

using namespace std;

#define OMEGA 0.3

int main(int argc, const char *argv[])
{
    ConfigManager configManager;
    Domain domain(2);
    Mesh mesh(domain, 2);
    Field<double, 2> u, v, f;
    Field<SpaceCoord> xd;
    TimeManager timeManager;
    IOManager io;
    int outputFileIdx;
    TimeLevelIndex<2> newIdx, halfIdx, oldIdx;

    Regrid regrid(mesh);

    double dt, dx, dy;
    int numIter;
    string outputPattern = "semi_lagrangian.%3s.nc";
    RegridMethod velocityInterpMethod;
    RegridMethod advectInterpMethod;

    if (argc != 2) {
        REPORT_ERROR("Configure file is needed!");
    }
    configManager.parse(argv[1]);
    dt = configManager.getValue("semi_lagrangian", "dt", 1.0);
    dx = configManager.getValue("semi_lagrangian", "dx", 0.01);
    dy = configManager.getValue("semi_lagrangian", "dy", 0.01);
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
    Time startTime(Date(2000, 1, 1), Seconds(0));
    Time endTime(Date(2000, 1, 1), Seconds(6*PI2/OMEGA));
    timeManager.init(startTime, endTime, dt);

    // Create velocity and density fields.
    u.create("u", "m s-1", "velocity component along x axis", mesh, CENTER, 2, true);
    v.create("v", "m s-1", "velocity component along y axis", mesh, CENTER, 2, true);
    f.create("f", "kg m-2", "tracer density", mesh, CENTER, 2);
    xd.create("xd", "m", "departure point coordinate", mesh, CENTER, 2);

    // Set initial guess for displacements.
    for (int j = mesh.js(FULL); j <= mesh.je(FULL); ++j) {
        for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
            xd(i, j).init(2);
            xd(i, j) = mesh.gridCoord(CENTER, i, j);
        }
    }

    // Set initial conditions.
    // NOTE: Velocity does not satisfy periodic boundary condition!
    SpaceCoord x0(2);
    x0.set(0.25, 0.5);
    newIdx = oldIdx+1;
    for (int j = mesh.js(FULL); j <= mesh.je(FULL); ++j) {
        for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
            const SpaceCoord &x = mesh.gridCoord(CENTER, i, j);
            double d = domain.calcDistance(x, x0);
            u(oldIdx, i, j) = -OMEGA*(x(1)-0.5);
            v(oldIdx, i, j) =  OMEGA*(x(0)-0.5);
            //u(oldIdx, i, j) = 0.01;
            //v(oldIdx, i, j) = 0.01;
            u(newIdx, i, j) = u(oldIdx, i, j);
            v(newIdx, i, j) = v(oldIdx, i, j);
            if (fabs(x(0)-x0(0)) >= 0.02 && d < 0.1) {
                f(oldIdx, i, j) = 1;
            } else if (fabs(x(0)-x0(0)) <= 0.02 && d < 0.1 && x(1) >= 0.55) {
                f(oldIdx, i, j) = 1;
            } else {
                f(oldIdx, i, j) = 0;
            }
        }
    }
    u.applyBndCond(oldIdx); u.applyBndCond(newIdx, true);
    v.applyBndCond(oldIdx); v.applyBndCond(newIdx, true);
    f.applyBndCond(oldIdx);

    // Set IO manager.
    io.init(timeManager);
    outputFileIdx = io.addOutputFile(mesh, outputPattern, Seconds(dt));
    io.addField(outputFileIdx, "double", FULL_DIMENSION, {&u, &v, &f});
    io.output<double, 2>(outputFileIdx, oldIdx, {&u, &v, &f});

    // Run the main loop.
    while (!timeManager.isFinished()) {
        halfIdx = oldIdx+0.5; newIdx = oldIdx+1;
        for (int j = mesh.js(FULL); j <= mesh.je(FULL); ++j) {
            for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
                const SpaceCoord &x0 = mesh.gridCoord(CENTER, i, j);
                SpaceCoord xm(2); // Middle point coordinate;
                MeshIndex meshIdx(2);
                Velocity U(2);
                // Calculate departure point.
                for (int iter = 0; iter < numIter; ++iter) {
                    xm() = x0()-0.5*(domain.diffCoord(x0, xd(i, j)));
                    domain.validateCoord(xm);
                    meshIdx.locate(mesh, xm);
                    regrid.run(velocityInterpMethod, halfIdx, u, xm, U(0), &meshIdx);
                    regrid.run(velocityInterpMethod, halfIdx, v, xm, U(1), &meshIdx);
                    meshIdx.locate(mesh, x0);
                    mesh.move(x0, -dt, U, meshIdx, xd(i, j));
                }
                // Calculate tracer mixing density on departure point (Note divergence is zero in this case!).
                meshIdx.locate(mesh, xd(i, j));
                regrid.run(advectInterpMethod, oldIdx, f, xd(i, j), f(newIdx, i, j), &meshIdx);
            }
        }
        f.applyBndCond(newIdx);
        timeManager.advance(); oldIdx.shift();
        io.output<double, 2>(outputFileIdx, oldIdx, {&u, &v, &f});
    }

    return 0;
}
