#include "geomtk/Cartesian.h"
#include <gsl/gsl_linalg.h>

using namespace std;

//#define USE_SUPERLU

int main(int argc, const char *argv[])
{
    ConfigManager configManager;
    Domain domain(1);
    Mesh mesh(domain);
    Field<double, 2> u, f;
    Field<double> fu;
#ifdef USE_SUPERLU
    arma::sp_mat A;
    arma::vec b, x;
#else
    gsl_vector *a1, *a2, *a3, *b, *x;
#endif
    TimeManager timeManager;
    IOManager io;
    int outputFileIdx;
    TimeLevelIndex<2> oldIdx, newIdx, halfIdx;

    double dt, dx;
    string outputPattern = "crank_nicolson.%3s.nc";

    if (argc != 2) {
        REPORT_ERROR("Configure file is needed!");
    }

    // Read configuration from file.
    configManager.parse(argv[1]);
    dt = configManager.getValue("crank_nicolson", "dt", 1.0);
    dx = configManager.getValue("crank_nicolson", "dx", 0.01);
    outputPattern = configManager.getValue("crank_nicolson", "output_pattern", outputPattern);

    // Set the one dimensional space axis.
    domain.setAxis(0, "x", "x axis", "m",
                   0, geomtk::BndType::PERIODIC,
                   1, geomtk::BndType::PERIODIC);

    // Set the discrete mesh on the domain.
    mesh.init(domain.axisSpan(0)/dx);

    // Set the time manager.
    Time startTime(Date(2000, 1, 1), Seconds(0));
    Time endTime(Date(2000, 1, 1), Seconds(1000));
    timeManager.init(startTime, endTime, dt);

    // Set up velocity and density fields.
    u.create("u", "m s-1", "velocity component along x axis", mesh, X_FACE, 1);
    f.create("f", "kg m-1", "tracer density", mesh, CENTER, 1);
    fu.create("fu", "kg s-1", "tracer mass flux", mesh, X_FACE, 1);

    // Create internal computing objects.
#ifdef USE_SUPERLU
    A.set_size(mesh.nx(FULL), mesh.nx(FULL));
    b.set_size(mesh.nx(FULL));
    x.set_size(mesh.nx(FULL));
#else
    a1 = gsl_vector_alloc(mesh.totalNumGrid(CENTER));
    a2 = gsl_vector_alloc(mesh.totalNumGrid(CENTER));
    a3 = gsl_vector_alloc(mesh.totalNumGrid(CENTER));
    b  = gsl_vector_alloc(mesh.totalNumGrid(CENTER));
    x  = gsl_vector_alloc(mesh.totalNumGrid(CENTER));
#endif

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
    outputFileIdx = io.addOutputFile(mesh, outputPattern, Seconds(dt));
    io.addField(outputFileIdx, "double", FULL_DIMENSION, {&f});
    io.output<double, 2>(outputFileIdx, oldIdx, {&f});

    // Run the main loop.
    double C = dt/dx*0.25;
    while (!timeManager.isFinished()) {
        newIdx = oldIdx+1; halfIdx = oldIdx+0.5;
        for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
            int j = i-1;
#ifdef USE_SUPERLU
            if (i == mesh.is(FULL)) {
                A(j, A.n_cols-1) = -C*u(newIdx, i-1);
                A(j, j+1       ) =  C*u(newIdx, i  );
            } else if (i == mesh.ie(FULL)) {
                A(j, j-1) = -C*u(newIdx, i-1);
                A(j, 0  ) =  C*u(newIdx, i  );
            } else {
                A(j, j-1) = -C*u(newIdx, i-1);
                A(j, j+1) =  C*u(newIdx, i  );
            }
            A(j, j) = 1+C*(u(newIdx, i)-u(newIdx, i-1));
            b(j) = f(oldIdx, i)-C*( u(oldIdx, i)                *f(oldIdx, i+1)+
                                   (u(oldIdx, i)-u(oldIdx, i-1))*f(oldIdx, i  )-
                                    u(oldIdx, i-1)*f(oldIdx, i-1));
#else
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
#endif
        }
#ifdef USE_SUPERLU
        x = spsolve(A, b);
#else
        gsl_linalg_solve_cyc_tridiag(a2, a3, a1, b, x);
#endif
        for (int i = mesh.is(FULL); i <= mesh.ie(FULL); ++i) {
            int j = i-1;
#ifdef USE_SUPERLU
            f(newIdx, i) = x(j);
#else
            f(newIdx, i) = gsl_vector_get(x, j);
#endif
        }
        f.applyBndCond(newIdx);
        timeManager.advance(); oldIdx.shift();
        io.output<double, 2>(outputFileIdx, oldIdx, {&f});
    }

#ifndef USE_SUPERLU
    gsl_vector_free(a1);
    gsl_vector_free(a2);
    gsl_vector_free(a3);
    gsl_vector_free(b);
    gsl_vector_free(x);
#endif

    return 0;
}
