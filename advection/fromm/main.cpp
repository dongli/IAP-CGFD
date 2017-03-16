#include "geomtk/Cartesian.h"

using namespace std;

int main(int argc, const char *argv[])
{
    ConfigManager configManager;
    Domain domain(1);
    Mesh mesh(domain, 2);
    Field<double, 2> u, f;
    Field<double> fu;
    TimeManager timeManager;
    IOManager io;
    int outputFileIdx;
    TimeLevelIndex<2> oldIdx, newIdx, halfIdx;

    double dt, dx, lw, bw;
    string outputPattern = "fromm.%3s.nc";

    if (argc != 2) {
        REPORT_ERROR("Configure file is needed!");
    }

    // Read configuration from file.
    configManager.parse(argv[1]);
    dt = configManager.getValue("fromm", "dt", 1.0);
    dx = configManager.getValue("fromm", "dx", 0.01);
    outputPattern = configManager.getValue("fromm", "output_pattern", outputPattern);

    // Set the one dimensional space axis.
    domain.setAxis(0, "x", "x axis", "m",
                   0, geomtk::BndType::PERIODIC,
                   1, geomtk::BndType::PERIODIC);

    // Set the discrete mesh on the domain.
    mesh.init(domain.axisSpan(0)/dx);

    // Set the time manager.
    Time startTime(Date(2000, 1, 1), Seconds(0));
    Time endTime(Date(2000, 1, 1), Seconds(200));
    timeManager.init(startTime, endTime, dt);

    // Set up velocity and density fields.
    u.create("u", "m s-1", "velocity component along x axis", mesh, X_FACE, 1, true);
    f.create("f", "kg m-1", "tracer density", mesh, CENTER, 1);
    fu.create("fu", "kg s-1", "tracer mass flux", mesh, X_FACE, 1);

    // Set the initial conditions.
    newIdx = oldIdx+1;
    for (int i = mesh.is(HALF); i <= mesh.ie(HALF); ++i) {
        u(oldIdx, i) = -0.005;
        u(newIdx, i) = -0.005;
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
    outputFileIdx = io.addOutputFile(mesh, outputPattern, Seconds(dt));
    io.addField(outputFileIdx, "double", FULL_DIMENSION, {&f});
    io.output<double, 2>(outputFileIdx, oldIdx, {&f});

    // Run the main loop.
    double C = dt/dx;
    while (!timeManager.isFinished()) {
        newIdx = oldIdx+1; halfIdx = oldIdx+0.5;
        for (int i = mesh.is(HALF); i <= mesh.ie(HALF); ++i) {
            // TODO: The following flux is only valid when u > 0.
            lw = 0.5*C*(      u(halfIdx, i)    *(f(oldIdx, i+1)+f(oldIdx, i))-
                        C*pow(u(halfIdx, i), 2)*(f(oldIdx, i+1)-f(oldIdx, i)));

            if (u(halfIdx, i) >= 0) {
              bw = 0.5*C*(      u(halfIdx, i)    *(3*f(oldIdx, i)-f(oldIdx, i-1))-
                          C*pow(u(halfIdx, i), 2)*(  f(oldIdx, i)-f(oldIdx, i-1)));

            } else {
              bw = 0.5*C*(      u(halfIdx, i)    *( -f(oldIdx, i+2)+3*f(oldIdx, i+1))-
                          C*pow(u(halfIdx, i), 2)*(  f(oldIdx, i+2)-  f(oldIdx, i+1)));
            }
            fu(i) = 0.5*(lw+bw);
        }
        fu.applyBndCond();
        for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
            f(newIdx, i) = f(oldIdx, i)-(fu(i)-fu(i-1));
        }
        f.applyBndCond(newIdx);
        timeManager.advance(); oldIdx.shift();
        io.output<double, 2>(outputFileIdx, oldIdx, {&f});
    }

    return 0;
}
