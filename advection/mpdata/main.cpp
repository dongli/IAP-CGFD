#include "geomtk/Cartesian.h"

#define DT 1
#define DX 0.01
#define IORD 4
#define EPS 1.0e-15
#define OUTPUT "mpdata.%3s.nc"

int main(int argc, const char *argv[])
{
    Domain domain(1);
    Mesh mesh(domain);
    Field<double, 2> u, f;
    Field<double> uc, fu, g;
    TimeManager timeManager;
    IOManager io;
    int outputFileIdx;
    TimeLevelIndex<2> oldIdx, newIdx, halfIdx;

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
    u.create("u", "m s-1", "velocity component along x axis", mesh, X_FACE, 1, true);
    f.create("f", "kg m-1", "tracer density", mesh, CENTER, 1);
    uc.create("uc", "m s-1", "antidiffusion velocity component along x axis", mesh, X_FACE, 1);
    fu.create("fu", "kg s-1", "tracer mass flux", mesh, X_FACE, 1);
    g.create("g", "kg m-1", "temporary tracer density", mesh, CENTER, 1);

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
        if (x(0) >= 0.05 && x(0) <= 0.25) {
            f(oldIdx, i) = 1.0;
        } else {
            f(oldIdx, i) = 0.0;
        }
    }
    f.applyBndCond(oldIdx);

    // Set up IO manager.
    io.init(timeManager);
    outputFileIdx = io.registerOutputFile(mesh, OUTPUT, geomtk::TimeStepUnit::STEP, 1);
    io.registerField(outputFileIdx, "double", FULL_DIMENSION, {&f});
    io.output<double, 2>(outputFileIdx, oldIdx, {&f});

    // Run the main loop.
    while (!timeManager.isFinished()) {
        newIdx = oldIdx+1; halfIdx = oldIdx+0.5;
        // Copy the initial velocity and tracer density.
        uc() = u(halfIdx); g() = f(oldIdx);
        for (int l = 1; l <= IORD; ++l) {
            if (l > 1) {
                // Calculate antidiffusion velocity.
                for (int i = mesh.is(HALF); i <= mesh.ie(HALF); ++i) {
                    uc(i) = (fabs(uc(i))*DX-pow(uc(i), 2)*DT)*
                            (g(i+1)-g(i))/(g(i+1)+g(i)+EPS)/DX;
                }
            }
            // Calculate the mass flux at cell interfaces.
            for (int i = mesh.is(HALF); i <= mesh.ie(HALF); ++i) {
                fu(i) = DT/DX*0.5*(     uc(i) *(g(i+1)+g(i))-
                                   fabs(uc(i))*(g(i+1)-g(i)));
            }
            fu.applyBndCond();
            for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
                g(i) = g(i)-(fu(i)-fu(i-1));
            }
            g.applyBndCond();
        }
        // Copy the final tracer density.
        f(newIdx) = g();
        timeManager.advance(); oldIdx.shift();
        io.output<double, 2>(outputFileIdx, oldIdx, {&f});
    }

    return 0;
}
