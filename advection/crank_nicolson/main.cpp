#include "geomtk/Cartesian.h"
#include <gsl/gsl_linalg.h>

#define DT 1
#define DX 0.01
#define OUTPUT "crank_nicolson.%3s.nc"

int main(int argc, const char *argv[])
{
    Domain domain(1);
    Mesh mesh(domain);
    Field<double, 2> u, f;
    Field<double> fu;
    gsl_vector *a1, *a2, *a3, *b, *x;
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
    u.create("u", "m s-1", "velocity component along x axis", mesh, X_FACE, 1);
    f.create("f", "kg m-1", "tracer density", mesh, CENTER, 1);
    fu.create("fu", "kg s-1", "tracer mass flux", mesh, X_FACE, 1);

    // Create internal computing objects.
    a1 = gsl_vector_alloc(mesh.totalNumGrid(CENTER));
    a2 = gsl_vector_alloc(mesh.totalNumGrid(CENTER));
    a3 = gsl_vector_alloc(mesh.totalNumGrid(CENTER));
    b  = gsl_vector_alloc(mesh.totalNumGrid(CENTER));
    x  = gsl_vector_alloc(mesh.totalNumGrid(CENTER));

    // Set the initial conditions.
    newIdx = oldIdx+1;
    for (int i = mesh.is(HALF); i <= mesh.ie(HALF); ++i) {
        u(oldIdx, i) = 0.005;
        u(newIdx, i) = 0.005;
    }
    u.applyBndCond(oldIdx);
    u.applyBndCond(newIdx);
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
    outputFileIdx = io.registerOutputFile(mesh, OUTPUT, geomtk::TimeStepUnit::STEP, 1);
    io.registerField(outputFileIdx, "double", FULL_DIMENSION, {&f});
    io.output<double, 2>(outputFileIdx, oldIdx, {&f});

    // Run the main loop.
    double C = DT/DX*0.25;
    while (!timeManager.isFinished()) {
        newIdx = oldIdx+1; halfIdx = oldIdx+0.5;
        for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
            int j = i-1;
            if (i == mesh.is(FULL)) {
                gsl_vector_set(a1, mesh.ie(FULL)-1, -C*u(newIdx, i-1));
            } else {
                gsl_vector_set(a1, j-1, -C*u(newIdx, i-1));
            }
            gsl_vector_set(a2, j, 1+C*(u(newIdx, i)-u(newIdx, i-1)));
            gsl_vector_set(a3, j, C*u(newIdx, i));
            gsl_vector_set(b, j, f(oldIdx, i)-C*( u(oldIdx, i)                *f(oldIdx, i+1)+
                                                 (u(oldIdx, i)-u(oldIdx, i-1))*f(oldIdx, i  )-
                                                               u(oldIdx, i-1) *f(oldIdx, i-1)));
        }
        gsl_linalg_solve_cyc_tridiag(a2, a3, a1, b, x);
        for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
            int j = i-1;
            f(newIdx, i) = gsl_vector_get(x, j);
        }
        f.applyBndCond(newIdx);
        timeManager.advance(); oldIdx.shift();
        io.output<double, 2>(outputFileIdx, oldIdx, {&f});
    }

    gsl_vector_free(a1);
    gsl_vector_free(a2);
    gsl_vector_free(a3);
    gsl_vector_free(b);
    gsl_vector_free(x);

    return 0;
}
