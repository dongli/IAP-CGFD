#include "geomtk/Cartesian.h"

using namespace std;

int main(int argc, const char *argv[])
{
    ConfigManager configManager;
    Domain domain(1);
    Mesh mesh(domain, 2);
    Field<double, 2> u, f;
    TimeManager timeManager;
    IOManager io;
    int outputFileIdx;
    TimeLevelIndex<2> oldIdx, newIdx, halfIdx;

    double dt, dx;
    string outputPattern = "beam_warming.%3s.nc";

    if (argc != 2) {
        REPORT_ERROR("Configure file is needed!");
    }

    // Read configuration from file.
    configManager.parse(argv[1]);
    dt = configManager.getValue("beam_warming", "dt", 1);
    dx = configManager.getValue("beam_warming", "dx", 0.01);
    outputPattern = configManager.getValue("beam_warming", "output_pattern", outputPattern);

    // Set the one dimensional space axis.
    domain.setAxis(0, "x", "x axis", "m",
                   0, geomtk::BndType::PERIODIC,
                   1, geomtk::BndType::PERIODIC);

    // Set the discrete mesh on the domain.
    mesh.init(domain.axisSpan(0)/dx);

    // Set the time manager.
    Time startTime(0*geomtk::TimeUnit::SECONDS);
    Time endTime(200*geomtk::TimeUnit::SECONDS);
    timeManager.init(startTime, endTime, dt);

    // Set up velocity and density fields.
    u.create("u", "m s-1", "velocity component along x axis", mesh, X_FACE, 1, true);
    f.create("f", "kg m-1", "tracer density", mesh, CENTER, 1);

    // Set the initial conditions.
    newIdx = oldIdx+1;
    for (int i = mesh.is(HALF); i <= mesh.ie(HALF); ++i) {
        u(oldIdx, i) = 0.005;
        u(newIdx, i) = 0.005;
    }
    u.applyBndCond(oldIdx);
    u.applyBndCond(newIdx, true);
    for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
        const SpaceCoord &x = mesh.gridCoord(CENTER, i);
        if (x(0) >= 0.05 && x(0) <= 0.1) {
            f(oldIdx, i) = 1.0;
        } else {
            f(oldIdx, i) = 0.0;
        }
    }
    f.applyBndCond(oldIdx);

    // Set up IO manager.
    io.init(timeManager);
    outputFileIdx = io.registerOutputFile(mesh, outputPattern, geomtk::TimeStepUnit::STEP, 1);
    io.registerField(outputFileIdx, "double", FULL_DIMENSION, {&f});
    io.output<double, 2>(outputFileIdx, oldIdx, {&f});

    // Run the main loop.
    while (!timeManager.isFinished()) {
        newIdx = oldIdx+1; halfIdx = oldIdx+0.5;
        for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
            double left = 0, right = 0;
            if (u(halfIdx, i-1) > 0) {
                double a = u(halfIdx, i-1)*dt/dx;
                left = dt/dx*u(halfIdx, i-1)*(0.5*(1-a)*f(oldIdx, i-2)-(2-a)*f(oldIdx, i-1)+0.5*(3-a)*f(oldIdx, i));
            }
            if (u(halfIdx, i) < 0) {
                double a = -u(halfIdx, i)*dt/dx;
                right = dt/dx*u(halfIdx, i)*(-0.5*(1+a)*f(oldIdx, i+2)+(2+a)*f(oldIdx, i+1)-0.5*(3+a)*f(oldIdx, i));
            }
            f(newIdx, i) = f(oldIdx, i)-(left+right);
        }
        f.applyBndCond(newIdx);
        timeManager.advance(); oldIdx.shift();
        io.output<double, 2>(outputFileIdx, oldIdx, {&f});
    }

    return 0;
}
