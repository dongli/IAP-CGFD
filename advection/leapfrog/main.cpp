#include "geomtk/Cartesian.h"

#define DT 1
#define DX 0.01
#define OUTPUT "leapfrog.%3s.nc"

int main(int argc, const char *argv[])
{
    Domain domain(1);
    Mesh mesh(domain);
    Field<double, 3> u, f;
    TimeManager timeManager;
    IOManager io;
    int outputFileIdx;
    TimeLevelIndex<3> l1, l2, l3;

    // Set the one dimensional space axis.
    domain.setAxis(0, "x", "x axis", "m",
                   0, geomtk::BndType::PERIODIC,
                   1, geomtk::BndType::PERIODIC);

    // Set the discrete mesh on the domain.
    mesh.init(domain.axisSpan(0)/DX);

    // Set the time manager.
    Time startTime(0*geomtk::TimeUnit::SECONDS);
    Time endTime(200*geomtk::TimeUnit::SECONDS);
    timeManager.init(startTime, endTime, DT);

    // Set up velocity and density fields.
    u.create("u", "m s-1", "velocity component along x axis", mesh, CENTER, 1);
    f.create("f", "kg m-1", "tracer density", mesh, CENTER, 1);

    // Set the initial conditions.
    l2 = l1+1; l3 = l2+1;
    for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
        u(l1, i) = 0.005;
    }
    u.applyBndCond(l1);
    u(l2) = u(l1); u(l3) = u(l1);
    for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
        const SpaceCoord &x = mesh.gridCoord(CENTER, i);
        if (x(0) >= 0.05 && x(0) <= 0.1) {
            f(l1, i) = 1.0;
        } else {
            f(l1, i) = 0.0;
        }
    }
    f.applyBndCond(l1);
    f(l2) = f(l1);

    // Set up IO manager.
    io.init(timeManager);
    outputFileIdx = io.registerOutputFile(mesh, OUTPUT, geomtk::TimeStepUnit::STEP, 1);
    io.registerField(outputFileIdx, "double", FULL_DIMENSION, {&f});
    io.output<double, 3>(outputFileIdx, l2, {&f});

    // Run the main loop.
    while (!timeManager.isFinished()) {
        l2 = l1+1; l3 = l2+1;
        for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
            f(l3, i) = f(l1, i)-DT/DX*(u(l2, i+1)*f(l2, i+1)-u(l2, i-1)*f(l2, i-1));
        }
        f.applyBndCond(l3);
        timeManager.advance(); l1.shift();
        io.output<double, 3>(outputFileIdx, l2, {&f});
    }

    return 0;
}
